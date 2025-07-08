sample="24WE0157_S30"

grep -v "*" ${sample}_gpu_HF.vcf  > ${sample}_gpu_HF_remove_star.vcf
docker run --rm -it -v $(pwd):/data tpesout/whatshap sh -c "whatshap phase -o ${sample}_gpu_phased.vcf --indels --no-reference ${sample}_gpu_HF_remove_star.vcf ${sample}_gpu.bam & chmod 777 ${sample}_gpu_phased.vcf"
