samples=$(cat sampleList.txt)
for sample in $samples;
do
	sudo nerdctl run --rm -v $(pwd):/workdir -v /purestorage/gene/Ref/hg19/:/refdir -w /workdir \
	geoffw/mosdepth:0.3.3 \
	sh -c "mosdepth --by /refdir/hg19_refseq_cds.bed --threads 8 --thresholds 1,10,20,50 -m -F 3844 --mapq 30 ${sample}/mosdepth ${sample}/nv_result/${sample}_gpu.bam &&
	zcat ${sample}/mosdepth.thresholds.bed.gz| \
	awk \"BEGIN{OFS=\\\"\\t\\\"} NR==1{print \\\"chrom\\\",\\\"start\\\",\\\"end\\\",\\\"region\\\",\\\"rate_1x\\\",\\\"rate_10x\\\",\\\"rate_20x\\\",\\\"rate_50x\\\";next} {len=\\\$3-\\\$2; print \\\$1,\\\$2,\\\$3,\\\$4,\\\$5/len,\\\$6/len,\\\$7/len,\\\$8/len}\"| \
	gzip -c > ${sample}/mosdepth.thresholds.ratio.bed.gz"

done
