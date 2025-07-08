files=$(cat $(pwd)/sampleList.txt)

for file_name in ${files}
do
        echo $file_name
        mkdir post_trim/
	
	docker run --rm -it -v $(pwd):/workdir staphb/trimmomatic \
	sh -c "java -jar /Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33  \
               /workdir/${file_name}_R1_001.fastq.gz \
               /workdir/${file_name}_R2_001.fastq.gz \
               /workdir/post_trim/${file_name}_R1_clean.fastq.gz \
               /workdir/post_trim/${file_name}_R1_unpaired.fastq.gz \
               /workdir/post_trim/${file_name}_R2_clean.fastq.gz \
               /workdir/post_trim/${file_name}_R2_unpaired.fastq.gz \
               ILLUMINACLIP:/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:2:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 -threads 32"

done

