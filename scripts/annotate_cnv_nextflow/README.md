## pull prerequisite docker images
* annotsv:
https://quay.io/repository/biocontainers/annotsv?tab=tags&tag=latest
---
## build docker through Dockerfile
```
docker build knotannotsv .
```
---
## run with nextflow pipeline
```
nextflow run run_annotation.nf --input_vcf [vcf_file] --output_dir [output_dir]
```
