########################################################################################################################
library(googlesheets)
library(data.table)
library(plyr)
library(dplyr)
library(stringr)
"%ni%" = Negate("%in%")

#######################################################################################################################

### Load cohort data
setwd('/ifs/res/taylorlab/chavans/gist')

msi<-fread('MSIscores.txt')
names(msi)[1] = "Tumor_Sample_Barcode"
names(msi)[4] = "MSIscore"

msi = msi %>% mutate(CMO_Sample_ID_fixed = str_replace_all(Tumor_Sample_Barcode,"s_","")) %>% 
  mutate(CMO_Sample_ID_fixed = str_replace_all(CMO_Sample_ID_fixed,"_","-")) %>%
  mutate(DMP_noIM = str_replace(CMO_Sample_ID_fixed,"-WES",""))

head(msi); head(msi); dim(msi)
length(unique(msi$Tumor_Sample_Barcode))

uniq_cmo_ids = unique(msi$DMP_noIM)

data_clinical_sample = fread("/ifs/res/taylorlab/chavans/dmp/mskimpact/data_clinical_sample.txt", skip = 4) %>% 
  filter(.,substr(SAMPLE_ID,1,13) %in% uniq_cmo_ids) %>% 
  select(PATIENT_ID,
         SAMPLE_ID,
         CANCER_TYPE,
         CANCER_TYPE_DETAILED,
         SAMPLE_TYPE,
         METASTATIC_SITE,
         PRIMARY_SITE,
         ONCOTREE_CODE) %>%
  mutate(Pool="Proj_06049_U")

data_clinical_patient = fread("/ifs/res/taylorlab/chavans/dmp/mskimpact/data_clinical_patient.txt", skip = 4) %>% 
  filter(.,substr(PATIENT_ID,1,9) %in% substr(uniq_cmo_ids,1,9)) %>% 
  select(PATIENT_ID,
         SEX,
         PARTC_CONSENTED_12_245) 
  #%>% mutate(Sex=str_replace_all(SEX,c("F"="Female","M"="Male"))) 

samplesheet=inner_join(data_clinical_patient,data_clinical_sample, by=c(PATIENT_ID = 'PATIENT_ID')) %>% 
  mutate(SAMPLE_ID_WES = str_c(substr(SAMPLE_ID,1,13),'-WES'))
head(samplesheet); dim(samplesheet) 

master=samplesheet 
dim(master) 
head(master)

########################################################################################################################
#FACETS CNA data 
cn_calls = fread('GIST_msk_exomes_hisens_cncf_GENELEVEL_070218.txt') %>% 
  mutate(DMP_noIM = str_replace_all(Tumor_Sample_Barcode,"s_","")) %>% 
  mutate(DMP_noIM = str_replace_all(DMP_noIM,"_","-")) %>%
  #mutate(DMP_noIM = str_replace_all(DMP_noIM,"-WES","")) %>%
  mutate(FACETS_CNA = ifelse(is.na(FACETS_CNA), NA, FACETS_CNA)) %>%
  dcast(Hugo_Symbol~DMP_noIM, value.var = 'FACETS_CNA')

head(cn_calls)
#length(unique(cn_calls$Tumor_Sample_Barcode))
length(names(cn_calls))

#names(cn_calls) = str_replace(names(cn_calls),'IM6','WES')
write.table(cn_calls, '/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/Proj_06049_U_data_CNA_n3.txt', sep="\t",quote=F,row.names=F)
########################################################################################################################
#FACETS WGD
wgd_facets<-fread('Proj_06049_U_WGD.txt') %>% 
  mutate(Sample_ID = str_replace_all(Sample_ID,"s_","")) %>% 
  mutate(Sample_ID = str_replace_all(Sample_ID,"_","-"))
head(wgd_facets)
unique(wgd_facets$Sample_ID) 

########################################################################################################################
#FACETS Segmentation

seg = fread('facets_log_ratio_.seg') %>% 
  mutate(V1 = str_replace_all(V1,"_","-"))
setnames(seg,c('ID','chrom','loc.start','loc.end','num.mark','seg.mean'))
dim(seg)
head(seg)
unique(seg$ID)
#seg = seg %>% mutate(ID = plyr::mapvalues(ID, master$CMO_Sample_ID_fixed, master$DMP_Sample_ID))
write.table(seg, '/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/Proj_06049_U_data_cna_hg19_n3.seg', sep="\t",quote=F,row.names=F)

########################################################################################################################
### Mutations
#cut -f1-109 
maf_ = fread('Proj_06049_U___SOMATIC.vep.filtered.facets.V3.postprocessed.mafanno.ccf.maf') %>%
  mutate(Tumor_Sample_Barcode = str_replace_all(Tumor_Sample_Barcode,"s_","")) %>% 
  mutate(Tumor_Sample_Barcode = str_replace_all(Tumor_Sample_Barcode,"_","-"))
  
exonic_mutations = c('Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins','Silent',
                     'Missense_Mutation','Nonsense_Mutation','Nonstop_Mutation','Splice_Site','Targeted_Region','Translation_Start_Site')

maf = filter(maf_, Variant_Classification %in% exonic_mutations)
dim(maf)
head(maf)
unique(maf$Tumor_Sample_Barcode) 

write.table(maf, '/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/Proj_06049_U_data_mutations_extended_n3.txt', sep="\t",quote=F,row.names=F)

########################################################################################################################

#Mut_Sig
sig = fread('signatures/Proj_06049_U_ListSignatures.txt') %>% 
  mutate(sample_name = str_replace_all(sample_name,"s_","")) %>% 
  mutate(sample_name = str_replace_all(sample_name,"_","-")) %>% 
  mutate(signatures=str_replace_all(signatures,",","|"),props=str_replace_all(props,",","|"),quasi_pvals=str_replace_all(quasi_pvals,",","|"))
head(sig)
unique(sig$sample_name) #rest did not have enough mutations

#Combine all
master_msi<-inner_join(master,msi,by=c(SAMPLE_ID_WES = 'CMO_Sample_ID_fixed'))
head(master_msi); dim(master_msi) #
master_msi_wgd<-inner_join(master_msi,wgd_facets,by=c(SAMPLE_ID_WES='Sample_ID'))
head(master_msi_wgd); dim(master_msi_wgd) #
master_msi_wgd_sig = left_join(master_msi_wgd,sig,by=c(SAMPLE_ID_WES='sample_name'))
head(master_msi_wgd_sig); dim(master_msi_wgd) #543

clin_sample = mutate(master_msi_wgd_sig) %>%
    select(PATIENT_ID, SAMPLE_ID = SAMPLE_ID_WES, CANCER_TYPE, CANCER_TYPE_DETAILED, SAMPLE_TYPE, METASTATIC_SITE, PRIMARY_SITE, ONCOTREE_CODE,
           FACETS_PURITY = Purity,
           FACETS_PLOIDY = Ploidy,
           GENOME_DOUBLED = WGD,
           FRACTION_CNA = FGA1,
           MSI_Score = MSIscore,
           Mutational_Signatures = signatures,
           Mutationa_Signatures_Proportions = props,
           Mutational_Signatures_Pval = quasi_pvals
           
    ) %>%
    mutate(PROJECT_CODE = 'Proj_06049_U (GIST)',
           INSTITUTE = 'MSKCC',
           CUD_CATEGORY = 'N/A',
           SAMPLE_CLASS = 'Tumor',
           LST = 'N/A',
           NTAI = 'N/A',
           HRD_LOH = 'N/A',
           BRCA_TYPE = 'N/A',
           BRCA_VARIANT = 'N/A',
           SOMATIC_BRCA_WT_STATUS = 'N/A'
    )
head(clin); dim(clin) #


clin_patient=mutate(master_msi_wgd_sig) %>%
           select(PATIENT_ID = PATIENT_ID,
                  SEX = SEX,
                  `12_245_PARTC_CONSENTED` = PARTC_CONSENTED_12_245)
head(clin_patient); dim(clin_patient)


write.table(clin, '/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/Proj_06049_U_data_clinical_sample_n3.txt',sep='\t',row.names=F,quote=F)
write.table(clin_patient, '/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/Proj_06049_U_data_clinical_patient_n3.txt',sep='\t',row.names=F,quote=F)

### Sample map

sample_map = master_msi_wgd_sig %>% 
  mutate(PATIENT_ID = PATIENT_ID, DMP_ID=SAMPLE_ID, SAMPLE_ID = SAMPLE_ID_WES) %>%
  select(PATIENT_ID, SAMPLE_ID, DMP_ID)
head(sample_map); dim(sample_map)
write.table(sample_map, '/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/Proj_06049_U_sample_map_n3.txt', sep="\t",quote=F,row.names=F)

##################################### END
