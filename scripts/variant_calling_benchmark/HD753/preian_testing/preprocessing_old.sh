main_folder=/home/cadilac/137_share/validation/TSO500
sampleList=$(cat sampleList.txt)
sampleList=HD753_8th
for i in ${sampleList}
do
	echo ${i}
	##remove CNV in vcf 
	bcftools view -i 'ALT!="<CNV>"' ${i}/main.vcf > ${i}/main.vcf_exclude_CNV.vcf
	bgzip ${i}/main.vcf_exclude_CNV.vcf
	bcftools index -t ${i}/main.vcf_exclude_CNV.vcf.gz
	##filter by target region
#	bcftools view ${i}/main.vcf_exclude_CNV.vcf.gz -R ${main_folder}/TST500C_manifest.bed > ${i}/main.vcf_exclude_CNV.vcf.gz_filter_by_menifest_Bed.vcf
#	bcftools view ${i}/main.vcf_exclude_CNV.vcf.gz -R ${main_folder}/TSO500_drop_unusual.bed > ${i}/main.vcf_exclude_CNV.vcf.gz_filter_by_drop_unusual_Bed.vcf
	bcftools view ${i}/main.vcf_exclude_CNV.vcf.gz -R ${main_folder}/TSO500_AGILENT_intersect.bed > ${i}/main.vcf_exclude_CNV.vcf.gz_filter_by_intersect_bed.vcf
	
	perl ~/137_share/147_backup/annovar/convert2annovar.pl -format vcf4 ${i}/main.vcf_exclude_CNV.vcf.gz_filter_by_intersect_bed.vcf -outfile ${i}/main.vcf_exclude_CNV.vcf.gz_filter_by_intersect_bed.vcf.avinput --includeinfo
done
