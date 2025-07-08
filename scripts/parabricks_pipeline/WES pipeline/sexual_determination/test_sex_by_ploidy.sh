## this script is used to determine a batch of samples' genders by the ploidies of chrX and chrY estimated by GATK germline CNV

for i in $(seq 0 35);
do
	#echo ${i};
	sampleID=$(cat contigPloidy/ploidy-calls/SAMPLE_${i}/sample_name.txt)
	ploidyX=$(cat contigPloidy/ploidy-calls/SAMPLE_${i}/contig_ploidy.tsv |grep chrX|cut -f2)
	ploidyY=$(cat contigPloidy/ploidy-calls/SAMPLE_${i}/contig_ploidy.tsv |grep chrY|cut -f2)
        #echo ${ploidyX} ${ploidyY}
        if [ ${ploidyX} -eq 2 ];then
		echo ${sampleID}"\tF\t"${ploidyX}"\t"${ploidyY}
	elif [ ${ploidyX} -eq 1 ] & [ ${ploidyY} -eq 1 ];then
		echo ${sampleID}"\tM\t"${ploidyX}"\t"${ploidyY}
	else
		echo ${sampleID}"\tfuck\t"${ploidyX}"\t"${ploidyY}
	fi
done >> test_sex_result.txt
