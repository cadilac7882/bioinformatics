ref_dir=/purestorage/gene/Ref/hg19
out_dir=$(pwd)
sampleList=$(cat ${out_dir}/sampleList.txt)
#sampleList=22W00692_S22

for file_name in ${sampleList}
do
        echo ${file_name}

        sudo nerdctl run --rm -v ${out_dir}:/outputdir \
                -v ${ref_dir}:/workdir \
                staphb/bcftools:1.21 \
                sh -c "bcftools +fill-tags /outputdir/${file_name}/nv_result/${file_name}_gpu.vcf -o  /outputdir/${file_name}/nv_result/${file_name}_gpu_bcfVAF.vcf -- -t VAF && \
		       bcftools norm -m-any /outputdir/${file_name}/nv_result/${file_name}_gpu_bcfVAF.vcf -o /outputdir/${file_name}/nv_result/${file_name}_gpu_bcfnorm.vcf && \
                       bcftools norm -f /workdir/ucsc.hg19.fasta /outputdir/${file_name}/nv_result/${file_name}_gpu_bcfnorm.vcf -o /outputdir/${file_name}/nv_result/${file_name}_gpu_bcfrealigned.vcf"


        sudo nerdctl run --rm -v ${out_dir}:/outputdir \
                -v ${ref_dir}:/workdir \
                broadinstitute/gatk:4.3.0.0 \
                sh -c "gatk SelectVariants -V /outputdir/${file_name}/nv_result/${file_name}_gpu_bcfrealigned.vcf -O /outputdir/${file_name}/nv_result/${file_name}_gpu_bcfrealigned_SNP.vcf -select-type SNP && \
                       gatk SelectVariants -V /outputdir/${file_name}/nv_result/${file_name}_gpu_bcfrealigned.vcf -O /outputdir/${file_name}/nv_result/${file_name}_gpu_bcfrealigned_INDEL.vcf -select-type INDEL && \
                       gatk VariantFiltration \
                                -V /outputdir/${file_name}/nv_result/${file_name}_gpu_bcfrealigned_SNP.vcf \
                                -filter \"QD < 2.0\" --filter-name \"QD2\" \
                                -filter \"QUAL < 30.0\" --filter-name \"QUAL30\" \
                                -filter \"SOR > 3.0\" --filter-name \"SOR3\" \
                                -filter \"FS > 60.0\" --filter-name \"FS60\" \
                                -filter \"MQ < 40.0\" --filter-name \"MQ40\" \
                                -filter \"MQRankSum < -12.5\" --filter-name \"MQRankSum-12.5\" \
                                -filter \"ReadPosRankSum < -8.0\" --filter-name \"ReadPosRankSum-8\" \
				--verbosity ERROR \
                                -O /outputdir/${file_name}/nv_result/${file_name}_gpu_bcfrealigned_SNP_filtered.vcf &&\
                        gatk VariantFiltration \
                                -V /outputdir/${file_name}/nv_result/${file_name}_gpu_bcfrealigned_INDEL.vcf \
                                -filter \"QD < 2.0\" --filter-name \"QD2\" \
                                -filter \"QUAL < 30.0\" --filter-name \"QUAL30\" \
                                -filter \"FS > 200.0\" --filter-name \"FS200\" \
                                -filter \"ReadPosRankSum < -20.0\" --filter-name \"ReadPosRankSum-20\" \
				--verbosity ERROR \
                                -O /outputdir/${file_name}/nv_result/${file_name}_gpu_bcfrealigned_INDEL_filtered.vcf &&\
                        gatk MergeVcfs -I /outputdir/${file_name}/nv_result/${file_name}_gpu_bcfrealigned_SNP_filtered.vcf \
                                -I /outputdir/${file_name}/nv_result/${file_name}_gpu_bcfrealigned_INDEL_filtered.vcf \
                                -O /outputdir/${file_name}/nv_result/${file_name}_gpu_HF.vcf \
				--VERBOSITY ERROR"

        rm -f ${file_name}/nv_result/${file_name}_gpu_bcf*

done

