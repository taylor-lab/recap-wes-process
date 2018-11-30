#!/usr/bin/env Rscript
#Usage Rscript ~/WES_QC_filters/scripts/run_WGD_test_all_batches.R ~/WES_QC_filters/facets/facets_mapping_.txt ~/WES_QC_filters/facets/Proj_93017_n581_WGD.txt

library(data.table)
source("WGD_function.R")

args<-commandArgs(trailingOnly = TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
    # default output file
    args[2] = "Summary_WGD.txt"
}

df<-args[1]
out<-args[2]

print(df)
print(out)

facets_mapping<-fread(df)
names(facets_mapping)[1]="Sample_ID"
names(facets_mapping)[2]="path"

#head(facets_mapping)
wgd_test(facets_mapping,out)

#Run mafanno
Rscript ~/git/facets-suite/mafAnno.R -m <filtered.maf> -f facets-mapping.txt -o <filtered.mafanno.maf>

#Run ccf (script included here)
Rscript ~/scripts/get_ccf_subclonality_all_batches.R <filtered.mafanno.maf> <filtered.mafanno.ccf.maf>
