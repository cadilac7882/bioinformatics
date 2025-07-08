#!/bin/bash

# Path to your BAM file
bam_file="$1"

# Path to your BED file
bed_file="$2"

echo $bam_file $bed_file
# Function to calculate mean MAPQ for a given region
calculate_mean_mapq() {
	local chrom="$1"
	local start="$2"
	local end="$3"
	#end=$(echo $end| tr -d '\r')
	local region=$(printf "%s:%d-%d" "$chrom" "$start" "$end")
	#echo "Debug: Calculating MAPQ for region $chrom:$start-$end"
	#echo "$bam_file"
	samtools view "$bam_file" "$chrom:$start-$end"| awk '{sum+=$5} END { if (NR > 0) meanMQ=sum/NR; else meanMQ=0; print "Mean MAPQ =",meanMQ",Reads =",NR }'
}

# Read regions from the BED file and calculate mean MAPQ for each
while IFS=$'\t' read -r chrom start end; do
	chrom=${chrom//$'\r'}
	start=${start//$'\r'}
	end=${end//$'\r'}
	
	result=$(calculate_mean_mapq ${chrom} ${start} ${end})


	echo "$chrom,$start,$end,$result"
done < "$bed_file"
