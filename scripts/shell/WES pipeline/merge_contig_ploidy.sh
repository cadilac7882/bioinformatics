#!/bin/bash

path=contigPloidy/ploidy-calls/

# Create the header line
header="CONTIG"
for i in $(seq 0 35); do
    sample_id=$(cat "${path}/SAMPLE_${i}/sample_name.txt")
    header+="\t${sample_id}"
done
echo -e "$header" > merged_contig_ploidy.tsv

# Process the first file to get the list of contigs
contigs=$(awk 'NR>2 {print $1}' "${path}/SAMPLE_0/contig_ploidy.tsv")

# Merge the data
while read contig; do
    echo -ne "$contig"
    for i in $(seq 0 35); do
        sample_file="${path}/SAMPLE_${i}/contig_ploidy.tsv"
        ploidy=$(awk -v cont="$contig" '$1 == cont && NR > 2 {print $2}' "$sample_file")
        echo -ne "\t${ploidy:-NA}"
    done
    echo ""
done <<< "$contigs" >> merged_contig_ploidy.tsv
