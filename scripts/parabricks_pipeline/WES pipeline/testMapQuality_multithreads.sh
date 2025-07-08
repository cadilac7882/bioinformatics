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
	nohup bash /scripts/testMapQuality.sh ${inputBam} ${resultDir}/target_p${j}.bed > ${resultDir}/tmp.mapq${j} &
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
	cat ${resultDir}/tmp.mapq${k} >> ${resultDir}/mybam.mapq
done
rm ${resultDir}/tmp.mapq* 
rm ${resultDir}/target_p*.bed

echo "complete!"
