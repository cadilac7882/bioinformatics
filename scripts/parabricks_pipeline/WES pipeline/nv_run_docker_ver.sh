ref_dir=/purestorage/gene/Ref/hg19
out_dir=$(pwd)
files=$(cat ${out_dir}/sampleList.txt)
#files=22W00692_S22

for file_name in ${files}
do
    echo ${file_name}
    output=${file_name}/nv_result
    mkdir ${out_dir}/${output}

	sudo nerdctl run --gpus '"device=10,11,12,13,14,15"' --rm -v ${ref_dir}:/workdir \
        -v ${out_dir}:/outputdir \
        -v ${out_dir}:/raid/myrun -w /workdir \
        nvcr.io/nvidia/clara/clara-parabricks:4.4.0-1 \
        pbrun germline --ref /workdir/ucsc.hg19.fasta \
              --in-fq /outputdir/${file_name}/${file_name}_R1_001.fastq.gz /outputdir/${file_name}/${file_name}_R2_001.fastq.gz \
              --knownSite /workdir/dbsnp_138.hg19.vcf.gz \
              --knownSite /workdir/1000G_phase1.indels.hg19.sites.vcf.gz \
              --knownSite /workdir/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz \
              --bwa-options="-M" \
			  --low-memory \
              --out-bam /outputdir/${output}/${file_name}_gpu.bam \
              --out-variants /outputdir/${output}/${file_name}_gpu.vcf \
              --out-recal-file /outputdir/${output}/${file_name}_gpu_recal.txt\
              --tmp-dir /outputdir
    rm ${file_name}/${file_name}_R*_001.fastq.gz 
done
