pb_dir=/purestorage/gene/Ref
out_dir=$(pwd)
files=$(cat ${out_dir}/sampleList.txt)
#files=22W00692_S22

for file_name in ${files}
do
        echo ${file_name}
        output=${file_name}/nv_result
        mkdir ${out_dir}/${output}

	sudo nerdctl run --gpus '"device=10,11,12,13,14,15"' --rm -v ${pb_dir}:/workdir \
		-v ${out_dir}:/outputdir \
		-v ${out_dir}:/raid/myrun -w /workdir \
		nvcr.io/nvidia/clara/clara-parabricks:4.4.0-1 \
		pbrun germline --ref /workdir/hg38/Homo_sapiens_assembly38.fasta \
                            --in-fq /outputdir/${file_name}/${file_name}_R1_001.fastq.gz /outputdir/${file_name}/${file_name}_R2_001.fastq.gz \
                            --knownSite /workdir/hg38/Homo_sapiens_assembly38.known_indels.vcf.gz \
                            --knownSite /workdir/hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
                            --bwa-options="-M" \
			    --read-group-sm ${file_name} \
                            --out-bam /outputdir/${output}/${file_name}_hg38.bam \
                            --out-variants /outputdir/${output}/${file_name}_hg38.vcf \
                            --out-recal-file /outputdir/${output}/${file_name}_hg38_recal.txt\
                            --tmp-dir /outputdir

done
