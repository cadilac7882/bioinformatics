sampleList=$(cat sampleList.txt)
#sampleList=22000466_S13
for sample in $sampleList
do	
	echo $sample
	gzip ${sample}/mybam.depth

done
