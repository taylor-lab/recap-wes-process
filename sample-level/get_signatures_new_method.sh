#!/bin/bash
#only if there are format issues, rerun vcf2maf for vep
#cmo_maf2maf --version 1.6.14 --vep-release 88 --input-maf data_mutations_unfiltered_93017_21batches.txt --outpput-maf data_mutations_unfiltered_93017_21batches.vep.maf

#--Edit this block
INPUT_MAF=Proj_06049_U___SOMATIC.vep.filtered.facets.V3.postprocessed.filter.maf
OUT_DIR=signatures
STUDY_ID="Proj_06049_U"
#-----------------

FILENAME=$(basename "$INPUT_MAF")
EXTENSION="${FILENAME##*.}"
PREFIX="${FILENAME%.*}"
echo $PREFIX
PDF_PLOT=${PREFIX}_mut_sig.pdf

#create snp-sorted input file

python /home/chavans/git/mutation-signatures/make_trinuc_maf.py \
  <(awk -F"\t" '{ if($1 == "Hugo_Symbol" || length($11) == 1 && length($13)==1)  print }' $INPUT_MAF) \
  $OUT_DIR/${PREFIX}_snpsorted.maf
 
#decompose
CI_CMD2="python /home/chavans/git/mutation-signatures/main.py \
  /home/chavans/git/mutation-signatures/Stratton_signatures30.txt \
  --spectrum_output $OUT_DIR/${PREFIX}_spectrum.txt \
  $OUT_DIR/${PREFIX}_snpsorted.maf \
  $OUT_DIR/${PREFIX}_snpsorted_decomposed.maf"

/common/lsf/9.1/linux2.6-glibc2.3-x86_64/bin/bsub \
-q sol -J mut_sig \
-cwd $OUT_DIR \
-e mut_sig1.stderr -o mut_sig1.stdout \
-We 00:59 -R "rusage[mem=5]" -M 10 -n 2 \
$CI_CMD2
  
#get significance  
CI_CMD3="~/git/signature.significance/signature.significance.sh --onepass -i $OUT_DIR/${PREFIX}_spectrum.txt -o $OUT_DIR/${PREFIX}_sigs.out"
/common/lsf/9.1/linux2.6-glibc2.3-x86_64/bin/bsub \
-q sol -J mut_sig \
-cwd $OUT_DIR \
-e mut_sig.stderr -o mut_sig.stdout \
-We 24:00 -R "rusage[mem=5]" -M 10 -n 20 \
$CI_CMD3