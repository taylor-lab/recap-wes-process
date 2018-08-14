library(data.table)
library(stringr)
library(reshape2)
library(dplyr)
"%nin%" <- Negate("%in%")

ngs_recurrent_fp = function(maf, hotspots = NULL){

if('TAG' %nin% names(maf))
   maf = maf %>% mutate(TAG = str_c(Chromosome,':',Start_Position,':',End_Position,':',Reference_Allele,':',Tumor_Seq_Allele2))

#Check for recurrenct variants prior to removing any variant // filtered in other sample, then remove
setnames(maf, "is-a-hotspot", "is_a_hotspot");
multi_tag = mutate(maf, patient = str_sub(Tumor_Sample_Barcode, 1, 10)) %>% 
                                          filter(.,is_a_hotspot != 'Y') %>% 
                                          distinct(patient, TAG, .keep_all = T) %>% 
                                          dplyr::count(TAG) %>% 
                                          filter(n>1)

filtered_vars = filter(maf, FILTER!='', FILTER!="PASS", FILTER!='RESCUE')

#Tagging ones in samples other than those where it's orignally tagged
recur_apply = filter(multi_tag, TAG %in% filtered_vars$TAG)
maf_flagged = maf %>% mutate(is_recurrent=ifelse(TAG %in% recur_apply$TAG,TRUE, FALSE))

print(dim(maf))
print(table(maf_flagged$is_recurrent))
print(dim(maf_flagged))

maf = maf_flagged
maf
}
