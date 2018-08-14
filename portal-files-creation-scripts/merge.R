#!/usr/bin/env Rscript
suppressMessages(library(data.table))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(kpjmisc))

########################################################################################################################
### Merge
### CUD
### Gray/BRCA
### ER BCa
### AKT1 exomes
########################################################################################################################

### CUD
cud_clin = fread('/ifs/work/bergerm1/camachon/P-06049/P-06049-Portal/data_clinical.txt') %>%
		filter(SAMPLE_ID != 'P-0001417-T01-WES')
cud_seg = fread('/ifs/work/bergerm1/camachon/P-06049/P-06049-Portal/CUD-FACETS-WES.seg') %>%
		filter(ID != 'P-0001417-T01-WES')
cud_maf = fread('/ifs/work/bergerm1/camachon/P-06049/P-06049-Portal/data_mutations_extended.txt') %>%
		filter(Tumor_Sample_Barcode != 'P-0001417-T01-WES')
cud_cna = fread('/ifs/work/bergerm1/camachon/P-06049/P-06049-Portal/data_CNA.txt') %>%
		dplyr::select(-`P-0001417-T01-WES`)

### BRCA/Gray
brca_clin = fread('/ifs/res/taylorlab/jonssonp/gray_brca/portal_instance/data_clinical.txt')
brca_seg = fread('/ifs/res/taylorlab/jonssonp/gray_brca/portal_instance/data_cna_hg19.seg')
brca_maf = fread('/ifs/res/taylorlab/jonssonp/gray_brca/portal_instance/data_mutations_extended.txt')
brca_cna = fread('/ifs/res/taylorlab/jonssonp/gray_brca/portal_instance/data_CNA.txt')

### ER BCa
erbca_clin = fread('/ifs/res/taylorlab/bandlamc/breast_mskcc/portal_instance/breast_endocrine/data_clinical.txt')
erbca_seg= fread('/ifs/res/taylorlab/bandlamc/breast_mskcc/portal_instance/breast_endocrine/data_cna_hg19.seg')
erbca_maf = fread('/ifs/res/taylorlab/bandlamc/breast_mskcc/portal_instance/breast_endocrine/data_mutations_extended.txt')
erbca_cna = fread('/ifs/res/taylorlab/bandlamc/breast_mskcc/portal_instance/breast_endocrine/data_CNA.txt')

### Ritika-given-files
prev_clin = fread('/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/previous_files/data_clinical_sample.txt', skip=4)
prev_clin_patient = fread('/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/previous_files/data_clinical_patient.txt', skip=4)
prev_seg = fread('/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/previous_files/mixed_impact_wes_2017_data_cna_hg19.seg')
prev_maf = fread('/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/previous_files/data_mutations_extended.txt')
prev_cna = fread('/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/previous_files/data_CNA.txt')

### Ovarian/bladder
proj6049r_clin = fread('/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/Proj_6049_R_data_clinical_sample.txt')
proj6049r_clin_patient = fread('/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/Proj_6049_R_data_clinical_patient.txt')
proj6049r_seg = fread('/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/Proj_06049_R_data_cna_hg19.seg')
proj6049r_maf = fread('/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/Proj_06049_R_data_mutations_extended.txt')
proj6049r_cna = fread('/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/Proj_06049_R_data_CNA.txt')

### Combine & Write
clin = rbind.fill(prev_clin, proj6049r_clin); head(clin); dim(clin)
clin_patient = rbind.fill(prev_clin_patient, proj6049r_clin_patient); head(clin_patient); dim(clin_patient)
seg = rbind(prev_seg, proj6049r_seg); head(seg); dim(seg)
cna = full_join(prev_cna, proj6049r_cna, by = 'Hugo_Symbol'); head(cna); dim(cna)
maf = rbind.fill(prev_maf, proj6049r_maf); head(maf); dim(maf)

write.table(clin, '/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/f.data_clinical_sample.txt',row.names=F, quote=F, sep="\t")
write.table(clin_patient, '/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/f.data_clinical_patient.txt',row.names=F, quote=F, sep="\t")
write.table(seg, '/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/f.data_cna_hg19.seg',row.names=F, quote=F, sep="\t")
write.table(maf, '/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/f.data_mutations_extended.txt',row.names=F, quote=F, sep="\t")
write.table(cna, '/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/f.data_CNA.txt',row.names=F, quote=F, sep="\t")

#Proj_93017
#proj93017_clin=fread('/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/Proj_93017_data_clinical_sample_010318.txt')
#proj93017_clin_ptn=fread('/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/Proj_93017_data_clinical_patient_010318.txt')
#proj93017_seg=fread('/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/Proj_93017_data_cna_hg19_010318.seg')
#proj93017_cna=fread('/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/Proj_93017_data_CNA_010318.txt')
#proj93017_maf=fread('/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/Proj_93017_data_mutations_extended_010318.txt')

#Proj_AKT
akt_clin=fread('/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/Proj_AKT_data_clinical_sample.txt')
akt_clin_ptn=fread('/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/Proj_AKT_data_clinical_patient.txt')
akt_seg=fread('/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/Proj_AKT_data_cna_hg19.seg')
akt_cna=fread('/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/Proj_AKT_data_CNA.txt')
akt_maf=fread('/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/Proj_AKT_data_mutations_extended.txt')

#So far
comb_prev_clin=fread('/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/all.data_clinical_sample_010318.txt')
comb_prev_clin_ptn=fread('/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/all.data_clinical_patient_010318.txt')
comb_prev_seg=fread('/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/all.data_cna_hg19_010318.seg')
comb_prev_cna=fread('/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/all.data_CNA_010318.txt')
comb_prev_maf=fread('/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/all.data_mutations_extended_010318.txt')

###NEW portal data #93017
proj93017_sample_map=fread('/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/Proj_93017_sample_map_n543_041018.txt') %>% mutate(DMP_WES=str_replace(DMP_ID,'IM6','WES')) %>% mutate(DMP_PAT=substr(DMP_ID, 1,9))
proj93017_clin=fread('/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/Proj_93017_data_clinical_sample_n543_041018.txt') %>% mutate(SAMPLE_ID=str_replace(SAMPLE_ID,'IM6','WES'))
proj93017_clin_ptn=fread('/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/Proj_93017_data_clinical_patient_n543_041018.txt') %>% mutate(SEX=ifelse(SEX=="F","Female","Male"))
proj93017_seg=fread('/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/Proj_93017_data_cna_hg19_n543_041018.seg') %>% mutate(ID=str_replace(ID,'IM6','WES'))
proj93017_cna=fread('/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/Proj_93017_data_CNA_n543_041018.txt')
proj93017_maf=fread('/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/Proj_93017_data_mutations_extended_n543_041018.txt')  %>% mutate(Tumor_Sample_Barcode=str_replace(Tumor_Sample_Barcode,'IM6','WES'))

####Ritika given updated portal files on 041018
assf_clin=fread('~/portal_instance_local/data_clinical_samples.txt') #788 24
assf_clin_ptn=fread('~/portal_instance_local/data_clinical_patient.txt') #744 3
assf_seg=fread('~/portal_instance_local/mixed_impact_wes_2017_data_cna_hg19.seg') #113683 6
assf_cna=fread('~/portal_instance_local/data_CNA.txt') #23782   734
assf_maf=fread('~/portal_instance_local/data_mutations_extended.txt') #354196    207

assf_clin_= assf_clin %>% filter(.,SAMPLE_ID %ni% proj93017_sample_map$DMP_WES)
assf_clin_ptn_ = assf_clin_ptn %>% filter(., PATIENT_ID %ni% proj93017_sample_map$DMP_PAT)
assf_seg_ = assf_seg %>% filter(., ID %ni% proj93017_sample_map$DMP_WES) #24399     6
assf_cna %>% select(-matches(paste0(proj93017_sample_map$DMP_WES, collapse = '|'))) -> assf_cna_
assf_maf_ = assf_maf %>% filter(.,Tumor_Sample_Barcode %ni% proj93017_sample_map$DMP_WES) #207860    207

### Combine & Write
clin = rbind.fill(proj93017_clin, assf_clin_); head(clin); dim(clin) #920
clin_patient = rbind.fill(proj93017_clin_ptn, assf_clin_ptn_); head(clin_patient); dim(clin_patient) #906
seg = rbind(proj93017_seg, assf_seg_); head(seg); dim(seg) # 71971     6
cna = full_join(proj93017_cna, assf_cna_, by = 'Hugo_Symbol'); head(cna); dim(cna) #23782   865
maf = rbind.fill(proj93017_maf, assf_maf_); head(maf); dim(maf); head(maf); dim(maf) #318049    249

###########

length(unique(clin$SAMPLE_ID))
unique(clin$SAMPLE_ID)
clin=mutate(clin, SAMPLE_ID=str_replace(SAMPLE_ID,substr(SAMPLE_ID,14,17),"-WES"))
unique(clin$SAMPLE_ID)

length(unique(clin_patient$PATIENT_ID))
unique(clin_patient$PATIENT_ID)

length(unique(seg$ID))
unique(seg$ID)
seg=mutate(seg, ID=str_replace(ID,substr(ID,14,17),"-WES"))
unique(seg$ID)

length(unique(names(cna)))
unique(names(cna))
names(cna)[2:length(names(cna))]=str_replace(names(cna)[2:length(names(cna))],substr(names(cna)[2:length(names(cna))],14,17),"-WES")
unique(names(cna))

length(unique(maf$Tumor_Sample_Barcode))
unique(maf$Tumor_Sample_Barcode)
maf=mutate(maf, Tumor_Sample_Barcode=str_replace(Tumor_Sample_Barcode,substr(Tumor_Sample_Barcode,14,17),"-WES"))
unique(maf$Tumor_Sample_Barcode)

write.table(clin, '/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/all.data_clinical_sample_akt_021418.txt',row.names=F, quote=F, sep="\t")
write.table(clin_patient, '/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/all.data_clinical_patient_akt_021418.txt',row.names=F, quote=F, sep="\t")
write.table(seg, '/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/all.data_cna_hg19_akt_021418.seg',row.names=F, quote=F, sep="\t")
write.table(cna, '/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/all.data_CNA_akt_021418.txt',row.names=F, quote=F, sep="\t")
write.table(maf, '/ifs/res/taylorlab/chavans/WES_QC_filters/portal_instance/all.data_mutations_extended_akt_021418.txt',row.names=F, quote=F, sep="\t")
