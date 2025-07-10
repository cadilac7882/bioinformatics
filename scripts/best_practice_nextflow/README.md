#--------------------------------------
# build docker through Dockerfile
#--------------------------------------
docker build -t dgm/secondary_analysis_toolkits:1.0 .

#--------------------------------------
# run with nextflow pipeline
#--------------------------------------
nextflow run best_practice_cpu.nf
