#!/usr/bin/env Rscript

#Usage : Rscript get_impact_maf.R
library(data.table)
library(plyr)
library(dplyr)
library(stringr)

  samplesheet<-fread("Proj_06049_U/r_001/Proj_06049_U_sample_pairing.txt") #DMP_request_93017_All_Batches.csv
  names(samplesheet) = c('CMO_N_ID','CMO_Sample_ID') 
  samplesheet<-samplesheet %>% mutate(CMO_Sample_ID=sub("s_","",CMO_Sample_ID))
  samplesheet = samplesheet %>% mutate(CMO_Sample_ID=gsub("_","-",CMO_Sample_ID), DMP = str_replace(CMO_Sample_ID,"-WES","")) 
  head(samplesheet)
  print(paste(length(unique(samplesheet$DMP))))

  #IMPACT maf : signedout #/ifs/res/taylorlab/chavans/data_mutations_extended.bak
      impact_maf_signed = fread('/home/chavans/dmp/mskimpact/data_mutations_extended.txt') %>%
      filter(substr(Tumor_Sample_Barcode,1,13) %in% samplesheet$DMP,Mutation_Status!='GERMLINE') %>%
      mutate(var_tag = str_c(Tumor_Sample_Barcode,':',Chromosome,':',Start_Position,':',End_Position,':',Reference_Allele,':',Tumor_Seq_Allele2), signed_out=TRUE)
      dim(impact_maf_signed)
      length(unique(impact_maf_signed$Tumor_Sample_Barcode))
  
  #Subset IMPACT maf for current cohort  
       write.table(distinct(impact_maf_signed),'GIST_impact.maf',sep="\t",row.names=FALSE,append=FALSE,quote=FALSE)

