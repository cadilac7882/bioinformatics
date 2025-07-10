## Build hap.py image through Dockerfile

```
## assume Dockerfile is downloaded in the current folder.
docker build -t dgm/hap.py .
```
---

## Compare somatic variant calling 

```
docker run --rm -v $(pwd):/data -w /data dgm/hap.py sh -c "som.py truth.vcf query.vcf -f confident.bed -o output_prefix -r reference.fa"
```


## Compare germline variant calling 
```
docker run --rm -v $(pwd):/data -w /data dgm/hap.py sh -c "hap.py truth.vcf query.vcf -f confident.bed -o output_prefix -r reference.fa"
```

## Benchmark with GIAB genome stratification for
[GIAB stratification ftp](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/)
```
docker run --rm -v $(pwd):/data -w /data dgm/hap.py sh -c "hap.py truth.vcf query.vcf -f confident.bed -o output_prefix -r reference.fa --stratification GRCh37@all/GRCh37-all-stratifications.tsv"
```

## Benchmark set for NA12878:
TN:
clinvar_20220320.vcf.gz	https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/archive_2.0/2022/
clinvar_20220320.vcf.gz.csi	index of clinvar_20220320.vcf.gz

TP:
HG001_GRCh37_1_22_v4.2.1_all.vcf.gz	https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh37/SupplementaryFiles/
HG001_GRCh37_1_22_v4.2.1_benchmark.bed	https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh37/
HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz	https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh37/
HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz.tbi	https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh37/
HG001_hg19_1_22_v4.2.1_benchmark.bed	contig with chr	altered from HG001_GRCh37_1_22_v4.2.1_benchmark.bed

## Benchmark set for NA12878:
HD75341130hg19genesdbsnpcosmic.vcf.gz download from Coa


### see https://github.com/Illumina/hap.py/tree/master for details.
