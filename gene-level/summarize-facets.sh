#!/bin/bash

#Usage: Run manually block by block 
### Get the preferred run for each sample
#ls s_C_*_d__s_C_*_d/facets_R0.5.6p1000c100/*purity.cncf.txt > purity-runs.txt

#while IFS=$'\t' read line
#do
#		SAMPLE=$(dirname $line)
#		RUN=$(basename $line)
#		REPLACE=$(grep $SAMPLE purity-runs.txt)
#		REPLACEMENT=$(sed "s/facets_c50p100/${RUN}/" <(echo $REPLACE))
#		sed -i "s,${REPLACE},${REPLACEMENT}," purity-runs.txt
#done < <(cat preferred-facets-run.txt)

PROJ_ID='GIST'

CNCF=($(cat purity-runs.txt |  tr "\n" "\t"))
OUT=($(cat purity-runs.txt | tr "\n" "\t" | sed 's/.cncf.txt/.out/g'))
SEG=($(cat purity-runs.txt | tr "\n" "\t" | sed 's/.cncf.txt/.seg/g'))

echo $CNCF
echo $OUT
echo $SEG

echo -e "Tumor_Sample_Barcode\tRdata_filename" > facets_mapping.txt
#paste <(cut -d"." -f8 purity-runs.txt | sed 's/printreads__//') <(sed 's/cncf.txt/Rdata/g' purity-runs.txt) >> facets_mapping.txt
#paste <(cut -d"/" -f9 purity-runs.txt | cut -d"_" -f1-5) <(sed 's/cncf.txt/Rdata/g' purity-runs.txt) >> facets_mapping.txt
paste <(cut -c1-19 purity-runs.txt) <(sed 's/cncf.txt/Rdata/g' purity-runs.txt) >> facets_mapping.txt

CNCF=($(echo ${CNCF[@]} | sed 's/purity/hisens/g'))
OUT=($(echo ${OUT[@]} | sed 's/purity/hisens/g'))
SEG=($(echo ${SEG[@]} | sed 's/purity/hisens/g'))

echo $CNCF
echo $OUT
echo $SEG

#check count of cncf files 
for index in ${!CNCF[@]}; do echo $index/${#CNCF[@]}; done
#check which files are missing.
for file in "${CNCF[@]}"; do     if [ ! -e "$file" ];     then echo "$file is missing";     fi; done

~/git/facets-suite/geneLevel.R -f ${CNCF[@]} -o ${PROJ_ID}_msk_exomes_hisens_cncf_GENELEVEL_070218.txt -m cncf

#MANUALLY fix this file to have only Tumor in the ID feild rather than both T-N
sed -i 's/__s_P_.*_WES//' GIST_msk_exomes_hisens_cncf_GENELEVEL_070218.txt
sed -i 's/_hisens//' GIST_msk_exomes_hisens_cncf_GENELEVEL_070218.txt

#Check uniq ids
cut -f1 GIST_msk_exomes_hisens_cncf_GENELEVEL_070218.txt | sort | uniq
