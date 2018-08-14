######purity-runs file
#ls /ifs/res/taylorlab/chavans/WES_QC_filters/res_93017_*/facets/s_C*_d*s_C_*_d/facets_R0.5.6s100n25c100p500/*purity.cncf.txt > purity-runs.txt
#ls /ifs/res/taylorlab/chavans/WES_QC_filters/res_93017_*/facets/s_C*_d*s_C_*_d*/facets_R0.5.6s100n25c100p500/*purity.cncf.txt >> purity-runs.txt

#!/bin/bash
OUT="purity-runs_.txt"
n=$(ls /ifs/res/taylorlab/chavans/WES_QC_filters/res_93017_*/facets/s_C*_d*s_C_*_d*/facets_R0.5.6s100n25c100p500/*purity.cncf.txt | wc -l)
start=$(ls /ifs/res/taylorlab/chavans/WES_QC_filters/res_93017_*/facets/s_C*_d*s_C_*_d*/facets_R0.5.6s100n25c100p500/*purity.cncf.txt | head)
ls /ifs/res/taylorlab/chavans/WES_QC_filters/res_93017_*/facets/s_C*_d*s_C_*_d*/facets_R0.5.6s100n25c100p500/*purity.cncf.txt > $OUT
echo "$OUT created with $n lines $start"

### Get the preferred run for each sample
#ls s_C_*_d__s_C_*_d/facets_R0.5.6p1000c100/*purity.cncf.txt > purity-runs.txt

while IFS=$'\t' read line
do
		SAMPLE=$(dirname $line);
		RUN=$(basename $line);
		echo $SAMPLE
		echo $RUN
		
		REPLACE=$(grep $SAMPLE purity-runs_.txt);
		echo $REPLACE
				
		REPLACEMENT=$(sed "s/facets_R0.5.6s100n25c100p500/${RUN}/" <(echo $REPLACE));
		echo $REPLACEMENT
		
		REPLACEMENT_DIR=$(dirname $REPLACEMENT);
		echo $REPLACEMENT_DIR

		REPLACEMENT2=$(sed "s/c500pc100_purity/_purity/" <(echo $REPLACEMENT));
		echo $REPLACEMENT2
		sed -i "s,${REPLACE},${REPLACEMENT2}," purity-runs_.txt;
done < <(cat preferred-facets-run.txt);

		  ##To add the suffix _purity to set of 5 files, add this bit to the loop above
		  #for file in ${REPLACEMENT_DIR}/*.*; :/1296
		  
		  #do V="$(echo $file |sed -e 's|_d\.|_d_purity\.|')"; 
  	  #echo $V; 
  		#mv $file $V; 
	    #done;

######log ratio file
#!/bin/bash
SEG=($(cat purity-runs_.txt | tr "\n" "\t" | sed 's/.cncf.txt/.seg/g'))
echo $SEG
cat ${SEG[@]} >> facets_log_ratio_.seg
sed -i '/^ID/d' facets_log_ratio_.seg
echo -e "ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean" > f.l.r_.seg
cat facets_log_ratio_.seg >> f.l.r_.seg

##Clean up the ID columns depending on the specific patterns of the sample IDs e.g. below, such that it only includes the Tumor Sample ID
sed -i 's/^s_\(P_.*_WES\)__.*_WES/\1/' facets_log_ratio_.seg
sed -i 's/_purity//' facets_log_ratio_.seg

#check the uniq id count 
#cut -f1 f.l.r_.seg | sort | uniq | wc -l

######facets estimates
#!/bin/bash
OUT_FILE='facets_estimates_.txt'
OUT=($(cat purity-runs_.txt | tr "\n" "\t" | sed 's/.cncf.txt/.out/g'))
echo $OUT
grep -E '(# Purity|# dipLogR)' ${OUT[@]}  > $OUT_FILE


#grep -E '(Purity|dipLogR)' ~/all_wes_batches/Proj_93017_*/r_001/facets/*purity.out > $OUT_FILE
