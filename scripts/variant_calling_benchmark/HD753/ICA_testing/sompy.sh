main_folder=/home/cadilac/137_share/validation/TSO500/HD753/
#sampleList=$(cat ${main_folder}/HD753/testing/sampleList.txt)
sampleList=VAL-02
for i in ${sampleList}
do
	mkdir ${main_folder}/HD753/testing/${i}/TP_menifest_som
	sudo docker run -it -v ${main_folder}:/data pkrusche/hap.py /opt/hap.py/bin/som.py /data/HD753/TP_hotspot/HD75341130hg19genesdbsnpcosmic.vcf.gz_filter_by_menifest_Bed.vcf /data/HD753/testing/${i}/main.vcf_exclude_CNV.vcf.gz_filter_by_menifest_Bed.vcf -f /data/TST500C_manifest.bed -o /data/HD753/testing/${i}/TP_menifest_som/result_test
	mkdir ${main_folder}/HD753/testing/${i}/TP_drop_unusual_som
	sudo docker run -it -v ${main_folder}:/data pkrusche/hap.py /opt/hap.py/bin/som.py /data/HD753/TP_hotspot/HD75341130hg19genesdbsnpcosmic.vcf.gz_filter_by_drop_unusual_Bed.vcf /data/HD753/testing/${i}/main.vcf_exclude_CNV.vcf.gz_filter_by_drop_unusual_Bed.vcf -f /data/TSO500_drop_unusual.bed -o /data/HD753/testing/${i}/TP_drop_unusual_som/result_test
done
