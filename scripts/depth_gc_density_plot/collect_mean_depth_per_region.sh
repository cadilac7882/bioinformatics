mapfile -t myArray < <(awk 'BEGIN{'\t'}{print $1" "$2}' sampleList.txt)

for row in "${myArray[@]}"
do
	sampleID=$(echo $row|cut -d' ' -f1|tr -d '\r')
	bamPath=$(echo $row|cut -d' ' -f2|tr -d '\r')
	echo ${sampleID}
	#echo ${bamPath}
	
	mkdir ${sampleID}
	ln -s ${bamPath} ${sampleID}/${sampleID}.bam

	time sudo nerdctl run --rm -it -u $(id -u):$(id -g) \
	-v /datalake_Intermediate/:/datalake_Intermediate \
	-v $(pwd):/workdir -v /purestorage/gene/Ref/:/Ref -w /workdir \
	broadinstitute/gatk:4.3.0.0 sh -c "samtools view -b -F 3332 -q 20 ${sampleID}/${sampleID}.bam --threads 8 > ${sampleID}/filtered.bam && \
   samtools index -b ${sampleID}/filtered.bam -@ 8 && \
   bedtools coverage -a /Ref/hg19/TruSeq_Exome_TargetedRegions_v1.2_hg19.bed -b ${sampleID}/filtered.bam -mean > ${sampleID}/bed_coverage_filtered.txt &&
   bedtools coverage -a /Ref/hg19/TruSeq_Exome_TargetedRegions_v1.2_hg19.bed -b ${sampleID}/${sampleID}.bam -mean > ${sampleID}/bed_coverage_original.txt &&
   bedtools coverage -a /workdir/TruSeq_Exome_TargetedRegions_v1.2_hg19.bed_binsize_150_window.bed_processed.bed -b ${sampleID}/${sampleID}.bam -mean > ${sampleID}/bed_coverage_original_bin150.txt &&
   bedtools coverage -a /workdir/TruSeq_Exome_TargetedRegions_v1.2_hg19.bed_binsize_150_window.bed_processed.bed -b ${sampleID}/filtered.bam -mean > ${sampleID}/bed_coverage_filtered_bin150.txt &&
   rm ${sampleID}/filtered.bam ${sampleID}/filtered.bam.bai"

	#time sudo nerdctl run --rm -it -v /datalake_Intermediate/:/datalake_Intermediate  -v $(pwd):/workdir -v /purestorage/gene/Ref/:/Ref -w /workdir \
        #broadinstitute/gatk:4.3.0.0 sh -c "bedtools coverage -a /workdir/TruSeq_Exome_TargetedRegions_v1.2_hg19.bed_binsize_150_window.bed_processed.bed -b ${sampleID}/${sampleID}.bam -mean > ${sampleID}/bed_coverage_original_bin150.txt && \
	#bedtools coverage -a /workdir/TruSeq_Exome_TargetedRegions_v1.2_hg19.bed_binsize_150_window.bed_processed.bed -b ${sampleID}/filtered.bam -mean > ${sampleID}/bed_coverage_filtered_bin150.txt"
done

