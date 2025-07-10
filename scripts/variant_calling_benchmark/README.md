## Build hap.py image through Dockerfile

```
## assume Dockerfile is downloaded in the current folder.
docker build -t dgm/hap.py .
```
---

## Compare somatic variant calling 


`docker run --rm -v $(pwd):/data -w /data dgm/hap.py sh -c "som.py truth.vcf query.vcf -f confident.bed -o output_prefix -r reference.fa"`


## Compare germline variant calling 

`docker run --rm -v $(pwd):/data -w /data dgm/hap.py sh -c "hap.py truth.vcf query.vcf -f confident.bed -o output_prefix -r reference.fa"`


### see https://github.com/Illumina/hap.py/tree/master for details.