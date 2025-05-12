#!/bin/bash

#============================================================#
# Preprocessing Script for Germline CNV Calling (WGS)
# Steps:
#   1. Generate intervals (no padding for WGS)
#   2. Download and prepare annotation tracks
#   3. Annotate intervals with GC, duplication, and mappability
# Requires: GATK 4.x, bedtools, awk, wget
#============================================================#

# Reference genome
REFERENCE="ucsc_hg19/ucsc.hg19.fasta"
FAI="${REFERENCE}.fai"

# Output names
TMP_INTERVALS="hg19.preprocessed.tmp.interval_list"
PROCESSED_INTERVALS="hg19.processed.interval_list"
ANNOTATED_OUTPUT="hg19.annotated.tsv"

#------------------------------------------------------------#
# Step 1: Generate intervals (padding = 0 for WGS)
#------------------------------------------------------------#
echo "Step 1: Preprocess intervals"

gatk PreprocessIntervals \
    -R $REFERENCE \
    --padding 0 \
    -imr OVERLAPPING_ONLY \
    -O $TMP_INTERVALS

# Filter intervals to retain only standard chromosomes (1–22, X, Y)
awk '$1 ~ /^@|chr([1-9]|1[0-9]|2[0-2]|X|Y)$/' $TMP_INTERVALS > $PROCESSED_INTERVALS

#------------------------------------------------------------#
# Step 2: Download and process annotation resources
#------------------------------------------------------------#

# 2.1 Segmental duplication BED file
SEGDUP_BED="hg19.nochr.SegDups.elements.no_gl.bed"

if [[ ! -f $SEGDUP_BED ]]; then
    echo "Downloading segmental duplication BED"
    wget -qO- https://storage.googleapis.com/gatk-best-practices/cnv_germline_pipeline/hg19.nochr.SegDups.elements.no_gl.bed | \
        awk '{print "chr"$0}' | \
        bedtools sort -faidx $FAI > $SEGDUP_BED
    gatk IndexFeatureFile -I $SEGDUP_BED
else
    echo "Segmental duplication BED already exists: $SEGDUP_BED"
fi

# 2.2 Mappability BED file
MAPPABILITY_BED="test_hg19.nochr.k100.umap.single.merged.bed"

if [[ ! -f $MAPPABILITY_BED ]]; then
    echo "Downloading mappability BED"
    wget -qO- https://storage.googleapis.com/gatk-best-practices/cnv_germline_pipeline/hg19.nochr.k100.umap.single.merged.bed.gz | \
        zcat | \
        awk '{print "chr"$0}' | \
        bedtools sort -faidx $FAI > $MAPPABILITY_BED
    gatk IndexFeatureFile -I $MAPPABILITY_BED
else
    echo "Mappability BED already exists: $MAPPABILITY_BED"
fi

#------------------------------------------------------------#
# Step 3: Annotate intervals with segmental duplication and mappability
#------------------------------------------------------------#
echo "Step 3: Annotating intervals"

gatk AnnotateIntervals \
    -L $PROCESSED_INTERVALS \
    -R $REFERENCE \
    --mappability-track $MAPPABILITY_BED \
    --segmental-duplication-track $SEGDUP_BED \
    -imr OVERLAPPING_ONLY \
    -O $ANNOTATED_OUTPUT

echo "Annotation complete: $ANNOTATED_OUTPUT"

# Additional resources: https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle
