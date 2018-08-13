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
	
  #NAME=${TBAM%%.*}
	
echo "$TBAM $NBAM";
  bsub -e $PWD -We 00:59 -n 10 -R "rusage[mem=10]" /home/jonssonp/git/msisensor/msisensor msi \
			-d /ifs/res/taylorlab/jonssonp/msk_glioma/impact/msisensor/impact/b37_dmp_microsatellites.list \
			-t $BAMDIR/'Proj_06049_U_indelRealigned_recal_'${TBAM}'.bam' \
			-n $BAMDIR/'Proj_06049_U_indelRealigned_recal_'${NBAM}'.bam' \
			-o ${TBAM}.msisensor\
  		-b 10
done < <(cat $BAMLIST)

#######
# Run manually to compile msi score from all samples

#grep -Ev "Total" *.msisensor | sed s/\.msisensor\:/"\t"/ > MSIscores.txt

