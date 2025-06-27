pb_dir=/purestorage/gene/parabricks/parabricks_sample/Ref/
script_path=/purestorage/gene/scripts/CNV/

sudo nerdctl run --rm -it -v $(pwd):/workdir -v ${script_path}:/script -v ${pb_dir}:/Ref geoffw/exomedepth1.1.10 \
	sh -c "Rscript /script/auto_exomeDepth_autosomal.R"
