inputBam=$1
targetBed=$2
devidedBy=$3
resultDir=$4

echo "Step1. Samtools depth."
nrows=$(cat $targetBed|wc -l)
echo "nrows=\t" $nrows
inter=$(awk -v nrows=$nrows -v devidedBy=$devidedBy 'BEGIN{printf("%d", nrows/devidedBy)}')
for i in `seq 1 $devidedBy`
do
        start=$(((i-1)*inter+1))
        if [ $i -ne $devidedBy ]
        then
                end=$((i*inter))
        else
                end=$nrows
        fi
        echo $start "-" $end
        sed -n "${start},${end}p" ${targetBed} > ${resultDir}/target_p${i}.bed
done


pids=""
for j in `seq 1 $devidedBy`
do
	nohup samtools depth -a -b ${resultDir}/target_p${j}.bed  ${inputBam} > ${resultDir}/tmp.depth${j} &
	#echo $i &
	pids="$pids $!"
done

for pid in $pids; do
	echo ${pid}
    ## cmd1 || cmd2 
    ## if cmd1 succeeds (i.e. exit code 0), then cmd2 will NOT be performed. 
    ## cmd2 will only be performed if cmd1 fails.
	wait ${pid}
	if [ $? -eq 0 ]; then
        	echo "SUCCESS - Job $pid exited with a status of $?"
 	else
		echo "FAILED - Job $pid exited with a status of $?"
	fi
done

for k in `seq 1 $devidedBy`
do
	cat ${resultDir}/tmp.depth${k} >> ${resultDir}/mybam.depth
done
rm ${resultDir}/tmp.depth* 
rm ${resultDir}/target_p*.bed

echo "Step2. Evaluate quality metrics."

result=${resultDir}/myreport.txt
# exclude 0x100 (secondary alignment)
total_PF_reads=$(samtools view  -F 0x100 -c ${inputBam} --threads 8)
echo "Total reads:\t" ${total_PF_reads} >> ${result}

duplicated_reads=$(samtools view -f 0x400 -c ${inputBam} --threads 8)
duplicated_rate=$(awk -v total_reads=$total_PF_reads -v duplicated_reads=$duplicated_reads 'BEGIN{printf("%.2f%%",duplicated_reads/total_reads*100)}')
echo "Duplicated rate:\t" ${duplicated_rate} >> ${result}

# exclude 0x100 (secondary alignment) and 0x4 (read unmapped)
mapped_reads=$(samtools view -F 0x100 -F 0x4 -c ${inputBam} --threads 8)
mapping_rate=$(awk -v total_reads=$total_PF_reads -v mapped_reads=$mapped_reads 'BEGIN{printf("%.2f%%",mapped_reads/total_reads*100)}')
echo "Mapping rate:\t" ${mapping_rate} >> ${result}

on_target_reads=$(samtools view -F 0x100 -F 0x4 -L ${targetBed} -c ${inputBam} --threads 8)
on_target_rate=$(awk -v total_reads=$total_PF_reads -v on_target_reads=$on_target_reads 'BEGIN{printf("%.2f%%",on_target_reads/total_reads*100)}')
echo "On target rate:\t" $on_target_rate >> ${result}

mean_depth=$(awk '{sum+=$3}END{printf("%.2f\n",sum/NR)}' ${resultDir}/mybam.depth)
echo "Mean depth:\t" $mean_depth >> ${result}

target_bases=$(awk '{sum+=($3-$2)}END{print sum}' ${targetBed})

uniformity=$(awk -v mean_depth=$mean_depth -v target_bases=$target_bases '{if($3>=mean_depth*0.2){cnt++}}END{printf("%.2f%%\n", cnt/target_bases*100)}' ${resultDir}/mybam.depth)
echo "Uniformity:\t" ${uniformity} >> ${result}

echo "complete!"
