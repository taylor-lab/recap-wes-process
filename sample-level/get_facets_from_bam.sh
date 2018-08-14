#!/bin/bash

### Run facets across all bams in mapped-bams.txt

BAMDIR=Proj_06049_U/r_001/alignments
BAMLIST=Proj_06049_U/r_001/Proj_06049_U_sample_pairing.txt

while IFS=$'\t' read line
do
  set $line
  TBAM=$2
  NBAM=$1
  echo "Tumor is $TBAM"
  #if [ $TBAM == 's_C_006117_P001_d' ] || [ $TBAM == 's_C_006118_P001_d' ] || [ $TBAM == 's_C_006119_P001_d' ]
  #then
  cmoflow_facets \
    --R_lib 0.5.6 \
    --normal-bam $BAMDIR/'Proj_06049_U_indelRealigned_recal_'${NBAM}'.bam' \
    --tumor-bam $BAMDIR/'Proj_06049_U_indelRealigned_recal_'${TBAM}'.bam' \
    --normal-name ${NBAM%.*} \
    --tumor-name ${TBAM%.*} \
    -c 100 \
    -pc 500
  #fi
  #--workflow-mode LSF \
    #echo "$TBAM $NBAM"
done < <(cat $BAMLIST)

