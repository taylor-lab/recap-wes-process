maf=/ifs/res/taylorlab/chavans/uterine/NEW/data_mutations_extended_108.vep.maf
oncotreecode="USARC"


  bsub -e $PWD -n 3 -We 0:59 \
    -R "rusage[mem=20] \
    -R 'select[internet]'" \
    "python ~/git/oncokb-annotator/MafAnnotator.py \
     -i $maf \
     -o ${maf}-oncokb-level.maf \
     -t $oncotreecode"
