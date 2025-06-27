mkdir vcf
out_dir=$(pwd)
files=$(cat ${out_dir}/sampleList.txt)

for file_name in ${files}
do
        #echo $file_name
	file_dir=${out_dir}/${file_name}
        cp ${out_dir}/${file_name}/nv_result/*_HF.vcf vcf/ 
done
