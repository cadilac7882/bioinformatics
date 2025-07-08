# Author: KD Yang
# Date: 2025-05-26
# Purpose: Generate PAR-masked chrY reference genome and related indices.
# Generally hg38 fasta and some hg19 fasta had masked the PAR of chrY, so we don't need to do this again.
# With masked PAR, all reads within PAR would be aligned to chrX without "secondary alignment", this will
# benefit the futher detection of SNV, CNV and SV.

# =======step 0.1: define bed file for Y-PAR ================ 
# hg19:
# chrY	10000	2649520
# chrY	59034049	59363566
# hg38:
# chrY	10000	2781479
# chrY	56887903	57217415
# Notice that bed file are 0-based
# =======step 0.2: check if PAR had been masked by N ========
sudo docker run --rm -it -v $(pwd):/refdir -w /refdir broadinstitute/gatk:4.3.0.0 \
	sh -c "samtools faidx ucsc.hg19.fasta chrY:10000-2649520"

# =======step 1: mask PAR through bedtools ==================
sudo docker run --rm -it -v $(pwd):/refdir -w /refdir broadinstitute/gatk:4.3.0.0 \
	sh -c "bedtools maskfasta -fi ucsc.hg19.fasta -bed hg19_Y_PAR_mask.bed -fo hg19_maskPAR/hg19_PAR_masked.fa"

# =======step 2.1: generate bwa indices =======================
sudo docker run --rm -it -v $(pwd):/refdir -w /refdir staphb/bwa \
	sh -c "bwa index hg19_maskPAR/hg19_PAR_masked.fa"

# =======step 2.2: generate gatk dictionary ===================
sudo docker run --rm -it -v $(pwd):/refdir -w /refdir broadinstitute/gatk:4.3.0.0 \
	sh -c "gatk CreateSequenceDictionary \
      R=hg19_maskPAR/hg19_PAR_masked.fa \
      O=hg19_maskPAR/hg19_PAR_masked.dict"

