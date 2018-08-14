#libraries & custom functions

suppressWarnings(library(data.table))
suppressWarnings(library(plyr))
suppressWarnings(library(dplyr))
suppressWarnings(library(stringr))
suppressWarnings(library(ggplot2))
suppressWarnings(library(reshape2))

specify_decimal = function(x, k) format(round(x, k), nsmall=k)
"%ni%" = Negate("%in%")

impact_maf_signed = fread('/ifs/res/taylorlab/chavans/WES_QC_filters/impact_data_mutations_extended_93017_n581.txt')

impact_maf_signed = mutate(impact_maf_signed, t_depth = t_alt_count + t_ref_count, t_var_freq = t_alt_count/t_depth)
unique(impact_maf_signed$Tumor_Sample_Barcode)

bam_dir = '/ifs/res/taylorlab/chavans/WES_QC_filters/all_bams'
bam_files=list.files(bam_dir,pattern="*.bam$",full.names=TRUE)
print(length(unique(bam_files)))


samplesheet<-fread("/ifs/res/taylorlab/chavans/WES_QC_filters/DMP_request_93017_581.csv") #DMP_request_93017_All_Batches.csv
samplesheet<-samplesheet %>% mutate(DMP=`DMP Sample ID`) %>% mutate(CMO_Sample_ID=gsub('-','_',`CMO Sample ID`)) 
samplesheet<-samplesheet %>% mutate(CMO_Sample_ID=sub("C_","s_C_",CMO_Sample_ID))
print(paste(length(unique(samplesheet$DMP))))

sample_pairing = fread('/ifs/res/taylorlab/chavans/WES_QC_filters/Sample_pairing_93017_n581.txt', header=FALSE)
names(sample_pairing)<-c("N","T") 
head(sample_pairing)#IMP verify that this is true for each batch#
sample_pairing<-filter(sample_pairing,N!="na" & T!="na")
print(sample_pairing[,2])
N=length(sample_pairing[,2])

maf_file='/ifs/res/taylorlab/chavans/WES_QC_filters/reprocess_variants/comb.processed.93017.581.recurr.postprocessed.postprocessed.filter.v2.maf'

exome_maf = fread(maf_file) %>% 
  filter(!grepl("N",Tumor_Sample_Barcode)) %>% 
  mutate(.,t_depth = t_alt_count + t_ref_count, t_var_freq = t_alt_count/t_depth)
print(dim(exome_maf))

#Add tags
exome_maf<-mutate(exome_maf, DMP = plyr::mapvalues(Tumor_Sample_Barcode, samplesheet$CMO_Sample_ID, samplesheet$DMP), 
                  var_tag = str_c(DMP,':',Chromosome,':',Start_Position,':',End_Position,':',Reference_Allele,':',Tumor_Seq_Allele2),
                  TAG = str_c(Chromosome,':',Start_Position,':',End_Position,':',Reference_Allele,':',Tumor_Seq_Allele2),
                  TAG_LC = str_c(Chromosome,':',Start_Position,':',Reference_Allele,':',Tumor_Seq_Allele2,':', Hugo_Symbol),
                  TAG_Gene = str_c(Chromosome, ':', Start_Position, ':', Hugo_Symbol))
print(dim(exome_maf))
print(unique(exome_maf$Tumor_Sample_Barcode))
print(unique(exome_maf$DMP))

impact_maf_signed = filter(impact_maf_signed, Tumor_Sample_Barcode %in% unique(samplesheet$DMP))
print(unique(impact_maf_signed$Tumor_Sample_Barcode))


if(length(intersect(unique(exome_maf$DMP),unique(impact_maf_signed$Tumor_Sample_Barcode)))==dim(samplesheet)[1]){
  print("no exome samples missing in impact maf")
  }

impact_maf_signed<-mutate(impact_maf_signed,
                          var_tag = str_c(substr(Tumor_Sample_Barcode,1,13),':',Chromosome,':',Start_Position,':',End_Position,':',Reference_Allele,':',Tumor_Seq_Allele2),
                          TAG = str_c(Chromosome,':',Start_Position,':',End_Position,':',Reference_Allele,':',Tumor_Seq_Allele2),
                          TAG_LC = str_c(Chromosome,':',Start_Position,':',Reference_Allele,':',Tumor_Seq_Allele2,':', Hugo_Symbol),
                          TAG_Gene = str_c(Chromosome, ':', Start_Position, ':', Hugo_Symbol))

N = dim(samplesheet)[1]

for(i in 1:N)
{ 
  DMP = as.character(samplesheet[i,'DMP']);  print(DMP)
  CMO = as.character(samplesheet[i,'CMO_Sample_ID']);  print(CMO)

  #subset maf to per sample
  impact_maf_signed_ = filter(impact_maf_signed,Tumor_Sample_Barcode==DMP)
  print(DMP)
  exome_maf_ = filter(exome_maf,Tumor_Sample_Barcode==CMO)
  print(CMO)

  #Compare mafs by var_tag and create exome-missed-variants-maf #unique(impact_maf_signed$var_tag) #unique(exome_maf$var_tag) 
  exome_missed_var_tags<-setdiff(unique(impact_maf_signed_$var_tag), unique(exome_maf_$var_tag))
  length(exome_missed_var_tags)
  exome_missed_impact_maf<-filter(impact_maf_signed_, var_tag %in% exome_missed_var_tags)
  print(paste("Done comparing IMPACT-maf to EXOME-maf..",length(unique(exome_missed_impact_maf$var_tag)),"count of exome misssed variants in",CMO, exome_missed_impact_maf$var_tag, exome_missed_impact_maf$TAG_Gene))

  #Run fillout
  if(length(unique(exome_missed_impact_maf$var_tag))>=1)
  {
    #Write out per sample mafs to fillout on (adds header as well)
    fillout_input_maf=paste0('/ifs/res/taylorlab/chavans/WES_QC_filters/reprocessed_variants_fillouts_v2/',CMO,'.per_sample_input_fo.maf')
    #print(CMO)
    #print(fillout_input_maf)
    write.table(exome_missed_impact_maf,fillout_input_maf,sep="\t",row.names=FALSE,append=FALSE,quote=FALSE)
    
    CMO_T = CMO
    CMO_N = plyr::mapvalues(CMO_T, sample_pairing$T, sample_pairing$N)
    tumor_bam=bam_files[as.numeric(grep(CMO_T,bam_files))]
    normal_bam=bam_files[as.numeric(grep(CMO_N,bam_files))]
    fillout_output_maf=paste0('/ifs/res/taylorlab/chavans/WES_QC_filters/reprocessed_variants_fillouts_v2/',CMO_T,'.per_sample_output_fo.maf')
    
    fillout_cmd=paste("cmo_fillout","-m",fillout_input_maf,"-b",tumor_bam,normal_bam,"-g GRCh37","-f 1","-p",fillout_output_maf,"-n 4","-v default")
    fillout_cmd_bsub=paste("bsub -e /ifs/res/taylorlab/chavans/WES_QC_filters/reprocessed_variants_fillouts_v2 -n 8 -R rusage[mem=5] -We 0:59",fillout_cmd)
    print(fillout_cmd_bsub)
    
    write.table(fillout_cmd_bsub,'/ifs/res/taylorlab/chavans/WES_QC_filters/reprocessed_variants_fillouts_v2/fillout_commands.txt',sep="\t",row.names=FALSE,append=TRUE,quote=FALSE)
    #system(fillout_cmd_bsub)
    #print(paste0("Running fill-out..",CMO_T))
  }#Run fillout
}#total_tumors

#Run the commands manually
#chmod 755 fillout_commands.txt
#.fillout_commands.txt

#Put together fillout mafs
fo_dir = '/ifs/res/taylorlab/chavans/WES_QC_filters/reprocessed_variants_fillouts_v2'
fo_files = list.files(fo_dir,pattern="*.per_sample_output_fo.maf$",full.names=TRUE)
print(length(unique(fo_files)))
N = length(unique(fo_files))
fillout_output_maf_fil_som = paste0(fo_dir,'/cohort_output_fo.T.maf')
                     
for(i in 1:N)
{
  #Filter fillout result file, and extract the read count data and remove the 'NORMAL' columns
                      curr_fo_file=fread(fo_files[i])
                      curr_fo_file=filter(curr_fo_file,!grepl("N",Tumor_Sample_Barcode)) #Consider removing only None
                      curr_fo_file_=curr_fo_file %>% mutate(n_depth=n_alt_count+n_ref_count, t_var_freq=t_alt_count/t_depth, fo_signed_out=TRUE)
                      head(curr_fo_file_)
                      dim(curr_fo_file_) 
                      write.table(distinct(curr_fo_file_),fillout_output_maf_fil_som,sep='\t',row.names=FALSE,quote=FALSE,append=TRUE)
                      print("Done filtering fill-out..")  
}

#Combine both exome and impact
cohort_fo = distinct(fread(fillout_output_maf_fil_som))
#DO THIS TO REMOVE MULTIPLE HEADERS -- perl -ne 'print unless $seen{$_}++' cohort_output_fo.T.maf > cohort_output_fo.T.uniq.maf
cohort_fo = fread('/ifs/res/taylorlab/chavans/WES_QC_filters/reprocessed_variants_fillouts_v2/cohort_output_fo.T.uniq.maf')
dim(cohort_fo)
exome_maf = exome_maf %>% mutate(fo_signed_out = FALSE)
comb_exome_impact = rbind.fill(exome_maf,cohort_fo)
write.table(comb_exome_impact,paste0(maf_file,'.fo.maf'),sep='\t',row.names=FALSE,quote=FALSE,append=FALSE)
print("Combined IMPACT-Exome file created..")

##combine coding only exome file with fo
coding_exome_maf = fread('/ifs/res/taylorlab/chavans/WES_QC_filters/reprocess_variants/comb.processed.93017.581.recurr.postprocessed.postprocessed.filter.v2.nona.coding.maf')
coding_exome_maf = coding_exome_maf %>% mutate(fo_signed_out = FALSE)
dim(coding_exome_maf)

comb_coding_exome_impact = rbind.fill(coding_exome_maf,cohort_fo)
dim(comb_coding_exome_impact)
write.table(comb_coding_exome_impact,'/ifs/res/taylorlab/chavans/WES_QC_filters/reprocess_variants/comb.processed.93017.581.recurr.postprocessed.filter.v2.nona.coding.fo.maf',sep='\t',row.names=FALSE,quote=FALSE,append=FALSE)

 