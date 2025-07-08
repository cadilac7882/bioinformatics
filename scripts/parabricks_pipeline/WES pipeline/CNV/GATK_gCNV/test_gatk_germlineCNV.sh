#!/bin/bash 
ref_dir=/purestorage/gene/Ref/hg19
sampleList=$(cat sampleList.txt) 

## step1: collect read counts from bam file
tmp_input_config=""
mkdir readCounts
for sample in ${sampleList}
do
	echo ${sample}	
	sudo nerdctl run --rm -v $(pwd):/outdir -v ${ref_dir}:/workdir broadinstitute/gatk:4.3.0.0 \
	sh -c "gatk CollectReadCounts -L /workdir/nextera-dna-targeted-regions-processed.interval_list \
						-R /workdir/ucsc.hg19.fasta \
						--interval-merging-rule OVERLAPPING_ONLY \
						-I /outdir/${sample}/nv_result/${sample}_gpu.bam \
						--format TSV \
						-O /outdir/readCounts/${sample}_counts.tsv"
						
	sed -i '27,94d' readCounts/${sample}_counts.tsv
	sed -i "s/sample/${sample}/g" readCounts/${sample}_counts.tsv
	
	tmp_input_config+=" -I /outdir/readCounts/"${sample}"_counts.tsv"
done					

## step2: filter problematic regions
sudo nerdctl run --rm -v $(pwd):/outdir -v ${ref_dir}:/workdir broadinstitute/gatk:4.3.0.0 \
	sh -c "gatk FilterIntervals -L /workdir/nextera-dna-targeted-regions-processed.interval_list \
	--annotated-intervals /workdir/nextera-dna-targeted-regions-processed.annotated.tsv \
	--exclude-intervals \"chr6:28477797-33448354\" \
	${tmp_input_config} \
	-imr OVERLAPPING_ONLY \
	-O /outdir/nextera-dna-targeted-regions-processed.cohort.GC.Mappability.SegDup.filtered.interval_list"

## step3: determine contig ploidy
mkdir contigPloidy
chmod 775 contigPloidy
echo ${tmp_input_config}
sudo nerdctl run --rm -v $(pwd):/outdir -v ${ref_dir}:/workdir broadinstitute/gatk:4.3.0.0 \
	sh -c "gatk DetermineGermlineContigPloidy -L /outdir/nextera-dna-targeted-regions-processed.cohort.GC.Mappability.SegDup.filtered.interval_list \
	--interval-merging-rule OVERLAPPING_ONLY \
	${tmp_input_config}  \
	--contig-ploidy-priors /workdir/gatk_GCNV_contig_ploidy_priors.tsv \
	--output /outdir/contigPloidy \
	--output-prefix ploidy \
	--verbosity DEBUG"

## step4: call CNV
mkdir gCNV_call
chmod 775 gCNV_call
sudo nerdctl run --rm -v $(pwd):/outdir -v ${ref_dir}:/workdir broadinstitute/gatk:4.3.0.0 \
	sh -c "gatk GermlineCNVCaller --run-mode COHORT \
			-L /outdir/nextera-dna-targeted-regions-processed.cohort.GC.Mappability.SegDup.filtered.interval_list \
			${tmp_input_config} \
			--annotated-intervals /workdir/nextera-dna-targeted-regions-processed.annotated.tsv \
			--contig-ploidy-calls /outdir/contigPloidy/ploidy-calls \
			--interval-merging-rule OVERLAPPING_ONLY \
			--output /outdir/gCNV_call \
			--output-prefix testGCNV \
			--verbosity DEBUG"


## step5: post-processing for CNV
mkdir gCNV_post
sampleIndex=$(ls gCNV_call/testGCNV-calls/|grep SAMPLE)
for sampleIndex in ${sampleIndex}
do
	index=$(echo ${sampleIndex}|cut -d'_' -f2)
	echo ${index}
	sampleID=$(cat gCNV_call/testGCNV-calls/${sampleIndex}/sample_name.txt)

	sudo nerdctl run --rm -v $(pwd):/outdir -v ${ref_dir}:/workdir broadinstitute/gatk:4.3.0.0 \
		sh -c "gatk PostprocessGermlineCNVCalls \
	        --model-shard-path /outdir/gCNV_call/testGCNV-model \
	        --calls-shard-path /outdir/gCNV_call/testGCNV-calls \
		--autosomal-ref-copy-number 2 \
	        --allosomal-contig chrX --allosomal-contig chrY \
	        --contig-ploidy-calls /outdir/contigPloidy/ploidy-calls \
		--sample-index ${index} \
		--output-genotyped-intervals /outdir/gCNV_post/${sampleID}_genotyped-intervals.vcf \
		--output-genotyped-segments /outdir/gCNV_post/${sampleID}_genotyped-segments.vcf \
		--output-denoised-copy-ratios /outdir/gCNV_post/${sampleID}_denoised_copy_ratios.tsv \
		--sequence-dictionary /workdir/ucsc.hg19.dict"
done
