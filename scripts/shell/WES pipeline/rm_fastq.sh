fastqList=$(find . -name "*.fastq.gz")
for fastq in ${fastqList}
do	
  echo "remove" ${fastq} !
  rm ${fastq}
done	
