## step1. make windows by bedtools makewindows
```
docker run --rm -v /home/bioinformatics/Reference/hg19/:/refdir \
 -v $(pwd):/workdir -w /workdir  dgm/secondary_analysis_toolkits:1.0 \
 sh -c "bedtools makewindows -b /refdir/TruSeq_Exome_TargetedRegions_v1.2_hg19.bed -w 150 -i srcwinnum > TruSeq_Exome_TargetedRegions_v1.2_hg19.bed_binsize_150_window.bed"
```
---
## step2. merge small windows 
```
./post_processing.sh TruSeq_Exome_TargetedRegions_v1.2_hg19.bed_binsize_150_window.bed TruSeq_Exome_TargetedRegions_v1.2_hg19.bed_binsize_150_window.processed.bed
```
---

## step3. annotate gc content by bedtools nuc
```
docker run --rm -v /home/bioinformatics/Reference/hg19/:/refdir \
 -v $(pwd):/workdir -w /workdir  dgm/secondary_analysis_toolkits:1.0 \
 sh -c "bedtools nuc -fi /refdir/ucsc.hg19.fasta -bed TruSeq_Exome_TargetedRegions_v1.2_hg19.bed_binsize_150_window.processed.bed > TruSeq_Exome_TargetedRegions_v1.2_hg19.bed_binsize_150_window.bed_processed.bed.gc_content.txt"
```
---
## step4. create sampleList.txt
```
foldername	path_of_bam 
25NGSA		/datalake_Intermediate/NextSeq2000/250326_VH00815_68_AACMMFLHV/25NGS-01_S29/nv_result/25NGS-01_S29_gpu.bam
```
---
## step5. create sampleList.txt
```
bash collect_mean_depth_per_region.sh
```
---

## step6. plot GC-depth density plot using R
```
Rscript generate_depth_gc_plot.R
```

