#!/bin/bash

# Function to print usage
print_usage() {
	echo "Usage: $0 [options]"
	echo "Options:"
	echo "  -b <file>    Input BAM file (optional)"
	echo "  -1 <file>    Input FASTQ R1 file (optional)"
	echo "  -2 <file>    Input FASTQ R2 file (optional)"
	echo "  -n <num>     Number of reads to select"
	echo "  -A <regA>    Targeted regions (e.g., chr1:100-200,chr2:300-400)"
	echo "  -B <regB>    Targeted regions (e.g., chr1:100-200,chr2:300-400)"
	echo "  -p <prefex>    Prefix of Output FASTQ file"
	echo "  -o <folder>    Output folder"
	exit 1
}

# Parse arguments
while getopts ":b:1:2:n:A:B:p:o:" opt; do
	case $opt in
		b) BAM_FILE="$OPTARG";;
		1) R1_FASTQ="$OPTARG";;
		2) R2_FASTQ="$OPTARG";;
		n) N_FRAGMENTS="$OPTARG";;
		A) REGION_A="$OPTARG";;
		B) REGION_B="$OPTARG";;
		p) OUTPUT_PREFIX="$OPTARG";;
		o) OUTPUT_FOLDER="$OPTARG";;
		\?) print_usage;;
	esac
done
if [ -z "$BAM_FILE" ] || [ -z "$N_FRAGMENTS" ] || [ -z "$OUTPUT_PREFIX" ] || [ -z "$REGION_A" ] || [ -z "$REGION_B" ] || [ -z "$R1_FASTQ" ] || [ -z "$R2_FASTQ" ] || [ -z "$OUTPUT_FOLDER" ]; then
	    print_usage
fi

if [ ! -d ${OUTPUT_FOLDER} ]; then
	mkdir ${OUTPUT_FOLDER}
	chmod 777 ${OUTPUT_FOLDER}
else
	echo "folder existed."
fi

chrA=$(echo $REGION_A|cut -d":" -f1)
chrB=$(echo $REGION_B|cut -d":" -f1)
docker run --rm -v $(pwd):/workdir -w /workdir staphb/samtools \
	sh -c "samtools view ${BAM_FILE} ${REGION_A}|awk -v chr=${chrB} 'BEGIN{FS=\"\t\"}{if(\$7==chr){print \$1}}'|sort|uniq > ${OUTPUT_FOLDER}/read_names_in_regionA.txt && \
	       samtools view ${BAM_FILE} ${REGION_B}|awk -v chr=${chrA} 'BEGIN{FS=\"\t\"}{if(\$7==chr){print \$1}}'|sort|uniq > ${OUTPUT_FOLDER}/read_names_in_regionB.txt && \
               cat ${OUTPUT_FOLDER}/read_names_in_regionA.txt ${OUTPUT_FOLDER}/read_names_in_regionB.txt |sort|uniq > ${OUTPUT_FOLDER}/all_read_names_in_region.txt && \
	       chmod 777 ${OUTPUT_FOLDER}/all_read_names_in_region.txt && \
	       rm ${OUTPUT_FOLDER}/read_names_in_regionA.txt ${OUTPUT_FOLDER}/read_names_in_regionB.txt && \
	       samtools view ${BAM_FILE} -N ${OUTPUT_FOLDER}/all_read_names_in_region.txt > ${OUTPUT_FOLDER}/reads_in_region.txt && \
	       chmod 777 ${OUTPUT_FOLDER}/reads_in_region.txt"

#shuf -n $N_READS reads_in_region.txt > selected_supporting_reads.txt
Rscript select_reads_v2.R ${N_FRAGMENTS} ${OUTPUT_FOLDER}/reads_in_region.txt ${OUTPUT_FOLDER}

#grep -F -x -v -f selected_supporting_reads.txt reads_in_region.txt|awk 'BEGIN{FS="\t"}{print $1}' > exclude_reads.txt

docker run --rm -v $(pwd):/workdir -w /workdir staphb/seqkit \
	sh -c "seqkit grep -f ${OUTPUT_FOLDER}/selected_reads.txt RNA_mt_a1_S2_R1_001.fastq.gz -o ${OUTPUT_FOLDER}/${OUTPUT_PREFIX}_sl_S2_R1_001.fastq.gz && \
	       seqkit grep -f ${OUTPUT_FOLDER}/selected_reads.txt RNA_mt_a1_S2_R2_001.fastq.gz -o ${OUTPUT_FOLDER}/${OUTPUT_PREFIX}_sl_S2_R2_001.fastq.gz && \
	       chmod 777 ${OUTPUT_FOLDER}/${OUTPUT_PREFIX}_sl_S2_R*_001.fastq.gz"

## -v: exclude, -f: select
#seqkit grep -v -f excluded_reads.txt RNA_mt_a1_S2_R1_001.fastq.gz -o ${OUTPUT_PREFIX}_S2_R1_001.fastq.gz
#seqkit grep -v -f excluded_reads.txt RNA_mt_a1_S2_R2_001.fastq.gz -o ${OUTPUT_PREFIX}_S2_R2_001.fastq.gz 
