main_folder=/home/cadilac/137_share/validation/TSO500
sampleList=$(cat sampleList)
#sampleList=VAL01
for i in ${sampleList}
do
	echo ${i}
	##remove CNV in vcf 
	#bcftools view -i 'FILTER=="PASS"' ${i}/${i}.hard-filtered.vcf > ${i}/${i}.hard-filtered-pass-only.vcf
	#bgzip ${i}/${i}.hard-filtered.vcf
	#bcftools index -t ${i}/${i}.hard-filtered.vcf.gz

	##filter by target region and pass-only
	#bcftools view ${i}/${i}.hard-filtered-pass-only.vcf.gz -R ${main_folder}/TSO500_AGILENT_intersect_notindifficult.bed > ${i}/${i}.hard-filtered-pass-only-TSO500_AGILENT_intersect_notindifficult.vcf
	#perl ~/137_share/147_backup/annovar/convert2annovar.pl -format vcf4 ${i}/${i}.hard-filtered-pass-only-TSO500_AGILENT_intersect_notindifficult.vcf -outfile ${i}/${i}.hard-filtered-pass-only-TSO500_AGILENT_intersect_notindifficult.vcf.avinput --includeinfo

	bcftools view ${i}/main.vcf_exclude_CNV.vcf.gz -R ${main_folder}/TSO500_AGILENT_intersect_notindifficult.bed > ${i}/main.vcf_exclude_CNV-TSO500_AGILENT_intersect_notindifficult.vcf
	perl ~/137_share/147_backup/annovar/convert2annovar.pl -format vcf4 ${i}/main.vcf_exclude_CNV-TSO500_AGILENT_intersect_notindifficult.vcf -outfile ${i}/main.vcf_exclude_CNV-TSO500_AGILENT_intersect_notindifficult.vcf.avinput --includeinfo


done
