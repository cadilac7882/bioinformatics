out_dir=$(pwd)
files=$(cat ${out_dir}/sampleList.txt)

for file_name in ${files}
do
        echo $file_name
	file_dir=${out_dir}/${file_name}
        mkdir ${file_dir}
	
        mv ${out_dir}/fastq/${file_name}_R*_001.fastq.gz ${file_dir}/

done
