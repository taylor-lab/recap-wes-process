#!/bin/bash


for maf in *_exome_impact_combined_raw.maf
do
#vepmaf=$line
  echo $maf
  NAME=$(basename $maf)
  extension="${NAME##*.}"
  filename="${NAME%.*}"
  sampleid=$(echo $maf | sed 's/_exome_impact_combined_raw.maf//')

  echo $filename
  echo $sampleid
#oncokb with cancer specific run for level 3B
     dmpid=$(grep $sampleid /home/chavans/WES_QC_filters/Sample_mapping_93017_n581.txt | cut -f1)
     oncotreecode=$(grep $dmpid /ifs/res/taylorlab/chavans/msk-impact/msk-impact/data_clinical.txt | cut -f15)
     bsub -e $PWD -n 3 -We 0:59 -R "rusage[mem=20] -R 'select[internet]'" "python ~/git/oncokb-annotator/MafAnnotator.py \
         -i $maf -o ${filename}-oncokb-level.maf -t $oncotreecode"
     #echo $dmpid
     #echo $oncotreecode
done

