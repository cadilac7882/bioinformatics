## prerequsite images and tools should be pulled
* ensemblorg/ensembl-vep:release_113.0  https://hub.docker.com/r/ensemblorg/ensembl-vep
* hydragenetics/automap:1.2  https://hub.docker.com/r/hydragenetics/automap
* Annovar  https://annovar.openbioinformatics.org/en/latest/user-guide/download/
---
## run with nextflow pipeline
```
nextflow run run_annotation.nf --input_vcf [vcf_file] --output_dir [output_dir]
```
