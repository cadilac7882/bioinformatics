#!/bin/bash
# Author: [Your Name]
# Date: 2025-05-08
# Purpose: Generate UCSC-style BED file (gene, exon, CDS) from NCBI RefSeq GTF
# Usage:
#   ./generate_refseq_bed.sh --build hg19 --feature CDS --output cds_hg19.bed [--select-only]

# ========== Argument Parsing ==========
BUILD=""
FEATURE=""
OUTPUT=""
SELECT_ONLY="false"
OUTDIR="output"

print_usage() {
    echo "Usage: $0 --build [hg19|hg38] --feature [gene|exon|CDS] --output output.bed [--select-only]"
    exit 1
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -b|--build)
            BUILD="$2"
            shift 2
            ;;
        -f|--feature)
            FEATURE="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT="$2"
            shift 2
            ;;
        -s|--select-only)
            SELECT_ONLY="true"
            shift
            ;;
        -h|--help)
            print_usage
            ;;
        *)
            echo "Unknown option: $1"
            print_usage
            ;;
    esac
done

# Validate arguments
if [[ -z "$BUILD" || -z "$FEATURE" || -z "$OUTPUT" ]]; then
    echo "Missing required arguments."
    print_usage
fi

mkdir -p "$OUTDIR"

# ========== Genome Build Setup ==========
if [[ "$BUILD" == "hg19" ]]; then
    GTF_URL="https://ftp.ncbi.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/105.20220307/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz"
    ACCESSION_TO_CHR_MAP=(
        'NC_000001.10=chr1' 'NC_000002.11=chr2' 'NC_000003.11=chr3' 'NC_000004.11=chr4'
        'NC_000005.9=chr5'  'NC_000006.11=chr6' 'NC_000007.13=chr7' 'NC_000008.10=chr8'
        'NC_000009.11=chr9' 'NC_000010.10=chr10' 'NC_000011.9=chr11' 'NC_000012.11=chr12'
        'NC_000013.10=chr13' 'NC_000014.8=chr14' 'NC_000015.9=chr15' 'NC_000016.9=chr16'
        'NC_000017.10=chr17' 'NC_000018.9=chr18' 'NC_000019.9=chr19' 'NC_000020.10=chr20'
        'NC_000021.8=chr21' 'NC_000022.10=chr22' 'NC_000023.10=chrX' 'NC_000024.9=chrY'
        'NC_012920.1=chrM'
    )
elif [[ "$BUILD" == "hg38" ]]; then
    GTF_URL="https://ftp.ncbi.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/110.20230912/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz"
    ACCESSION_TO_CHR_MAP=(
        'NC_000001.11=chr1' 'NC_000002.12=chr2' 'NC_000003.12=chr3' 'NC_000004.12=chr4'
        'NC_000005.10=chr5' 'NC_000006.12=chr6' 'NC_000007.14=chr7' 'NC_000008.11=chr8'
        'NC_000009.12=chr9' 'NC_000010.11=chr10' 'NC_000011.10=chr11' 'NC_000012.12=chr12'
        'NC_000013.11=chr13' 'NC_000014.9=chr14' 'NC_000015.10=chr15' 'NC_000016.10=chr16'
        'NC_000017.11=chr17' 'NC_000018.10=chr18' 'NC_000019.10=chr19' 'NC_000020.11=chr20'
        'NC_000021.9=chr21' 'NC_000022.11=chr22' 'NC_000023.11=chrX' 'NC_000024.10=chrY'
        'NC_012920.1=chrM'
    )
else
    echo "Unsupported genome build: $BUILD"
    exit 1
fi

# Validate feature type
if [[ "$FEATURE" != "gene" && "$FEATURE" != "exon" && "$FEATURE" != "CDS" ]]; then
    echo "Unsupported feature type: $FEATURE"
    echo "Must be one of: gene, exon, CDS"
    exit 1
fi

# ========== Step 1: Download GTF ==========
GTF_FILE="$OUTDIR/${GTF_URL##*/}"
if [[ ! -f "$GTF_FILE" ]]; then
    echo "Downloading GTF for $BUILD..."
    wget -O "$GTF_FILE" "$GTF_URL"
else
    echo "GTF already exists: $GTF_FILE"
fi

# ========== Step 2: Extract Features to BED ==========
BED_TMP="$OUTDIR/temp_${FEATURE}.bed"
echo "Extracting feature: $FEATURE (RefSeq Select = $SELECT_ONLY)..."

zcat "$GTF_FILE" \
| awk -v feature="$FEATURE" -v select="$SELECT_ONLY" '$3 == feature {
    if (select == "true" && $0 !~ /tag "RefSeq Select"/) next;
    print $0;
}' \
| awk -F '\t' -v feat="$FEATURE" 'BEGIN{OFS="\t"} {
    gene="."; transcript="."; exon=".";
    split($9, attrs, ";");
    for (i in attrs) {
        gsub(/^ +| +$/, "", attrs[i]);
        if (attrs[i] ~ /^gene /) {
            match(attrs[i], /"[^"]+"/);
            gene = substr(attrs[i], RSTART+1, RLENGTH-2);
        }
        if (attrs[i] ~ /^transcript_id /) {
            match(attrs[i], /"[^"]+"/);
            transcript = substr(attrs[i], RSTART+1, RLENGTH-2);
        }
        if (attrs[i] ~ /^exon_number /) {
            match(attrs[i], /"[^"]+"/);
            exon = substr(attrs[i], RSTART+1, RLENGTH-2);
        }
    }
    label = (feat=="gene" ? gene :
             feat=="exon" ? gene "_" transcript "_exon" exon :
             feat=="CDS"  ? gene "_" transcript "_cds_exon" exon : ".");
    print $1, $4-1, $5, label, ".", $7;
}' > "$BED_TMP"

# ========== Step 3: Convert NC_ to chrN ==========
echo "Mapping RefSeq accessions to UCSC chr names..."
cp "$BED_TMP" "$OUTPUT"

for map in "${ACCESSION_TO_CHR_MAP[@]}"; do
    from=${map%%=*}
    to=${map##*=}
    sed -i "s/^$from/$to/" "$OUTPUT"
done

echo "BED file generated: $OUTPUT"
