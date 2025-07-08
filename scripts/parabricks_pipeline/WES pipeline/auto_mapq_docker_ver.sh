ref_dir=/purestorage/gene/Ref/hg19
script_dir=$(pwd)/scripts
sampleList=$(cat sampleList.txt)
#sampleList=22W00690T_S12
for sample in ${sampleList}
do	
	echo $sample
	#mkdir ${sample}
	start=`date +%s.%3N`

	sudo nerdctl run --rm -it -v $(pwd):/workdir -v ${ref_dir}:/ref_dir -v ${script_dir}:/scripts staphb/samtools:1.21 \
	sh -c "sh /scripts/testMapQuality_multithreads.sh /workdir/${sample}/nv_result/${sample}_gpu.bam  /ref_dir/nextera-dna-targeted-regions.bed 12 /workdir/${sample}"

	end=`date +%s.%3N`
	runtime=$( echo "$end - $start" |bc)
	hours=$( echo "$runtime / 3600" |bc)
	minutes=$( echo "($runtime % 3600) / 60" |bc)
	seconds=$( echo "($runtime % 3600) % 60" |bc)
	LC_NUMERIC=C printf "Runtime: %dh:%dm:%.3fs\n" $hours $minutes $seconds

done

