ref_dir=/purestorage/gene/Ref/hg19
script_dir=/purestorage/gene/scripts
sampleList=$(cat sampleList.txt)
#sampleList=22W00690T_S12
for sample in ${sampleList}
do	
	echo $sample
	#mkdir ${sample}
	start=`date +%s.%3N`

	sudo nerdctl run --rm -it -v $(pwd):/workdir -v ${ref_dir}:/ref_dir -v ${script_dir}:/scripts staphb/samtools:1.21 \
	sh -c "sh /scripts/samtools_metric_multithreads.sh /workdir/${sample}/nv_result/${sample}_gpu.bam  /ref_dir/nextera-dna-targeted-regions.bed 12 /workdir/${sample}"

	end=`date +%s.%3N`
	runtime=$( echo "$end - $start" |bc)
	hours=$( echo "$runtime / 3600" |bc)
	minutes=$( echo "($runtime % 3600) / 60" |bc)
	seconds=$( echo "($runtime % 3600) % 60" |bc)
	LC_NUMERIC=C printf "Runtime: %dh:%dm:%.3fs\n" $hours $minutes $seconds

done


touch report_summary.csv
echo "ID,Total reads,Duplicated rate,Mapping rate,On target rate,Mean depth,Uniformity,QC" >> report_summary.csv

for sample in ${sampleList}
do
        echo ${sample}
        sed 's/%//' ${sample}/myreport.txt| \
        awk -v var=${sample} 'BEGIN{FS="\t "}{ for(i=1; i<=NF; i++){ a[NR,i] = $i }} \
             NF>p { p = NF;} \
             END { for(j=2; j<=p; j++) { str=a[1,j]; for(i=2; i<=NR; i++){if(i==5){str=str","a[i,j]}else{str=str","a[i,j]"%"};} \
                  if(a[1,2]>=30000000 && a[3,2]>=95 && a[4,2]>=40 && a[5,2]>=50 && a[6,2]>=90){q="PASS"}else{q="FAIL"} \
                  print var","str","q \
               }}' >> report_summary.csv
done

