scriptsDir=/purestorage/gene/scripts

##
echo $(pwd)

## 1. generate sampleList
echo "1. generate sampleList:"

ls fastq|grep "fastq.gz"|awk 'BEGIN{FS="_R"}{print $1}'|uniq|grep -v "Undetermine" > sampleList.txt

sampleCnt=$(cat sampleList.txt|wc -l)
if [ "$sampleCnt" -lt 2 ];
then
    echo "${sampleCnt} sample was detected"
else
    echo "${sampleCnt} samples were detected"
fi
echo "------------------------------------"

## 2. archive fastqs
echo "2. archive fastq by sample:"
sh archive_fastq.sh

echo "------------------------------------"
## 3. run secondary analysis by parabricks
echo "3. secondary analysis by parabricks:"
sh nv_run_docker_ver.sh

echo "------------------------------------"
## 4. run hard filtering by gatk
echo "4. hard filtering by gatk:"
sh hardfiltering_docker_ver.sh 

echo "------------------------------------"
## 5. generate QC metrics
echo "5. generate QC report:"
sh auto_QC_metric_docker_ver.sh

echo "------------------------------------"
## 6. mileup compression
echo "6. mileup compression:"
sh auto_zip.sh
echo "------------------------------------"

## 7. call CNV by GATK germlineCNV caller (cohort mode)
bash test_gatk_germlineCNV.sh
bash merge_contig_ploidy.sh

## 8. mosdepth
bash mosdepth_cds.sh

