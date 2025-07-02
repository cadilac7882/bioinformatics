#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ========================================================================================
// PARAMETER DEFINITIONS
// ========================================================================================
params.input_vcf = false // Input VCF file. Mandatory parameter.
params.output_dir = false // Default output directory

if (!params.input_vcf || !params.output_dir) {
    exit 1, "ERROR: Both --input_vcf and --output_dir must be specified."
}

log.info """
         V C F   A N N O T A T I O N   P I P E L I N E
         ================================================
         Input VCF     : ${params.input_vcf}
         Output Dir    : ${params.output_dir}
         """
// ========================================================================================
// WORKFLOW DEFINITION
// ========================================================================================

workflow {
   Channel
        .fromPath(params.input_vcf)
        .map { vcf_file ->
            def id = vcf_file.simpleName // e.g., "sample1" from "sample1.vcf"
            def outdir = file(params.output_dir) // Convert output dir string to a path object
            return [ id, vcf_file, outdir ]
        }
        .view()
        .set { ch_input_vcf }

    // 1. run VEP annotation
    VEP_ANNOTATION(ch_input_vcf)
    VEP_ANNOTATION.out.vep_json_ch.view { id, json ->
        "SUCCESS (VEP): annotation for sample '$id' is complete. Output at: $json"
    }
    
	// 2. run annovar annotation
    ANNOVAR_ANNOTATION(ch_input_vcf)
    ANNOVAR_ANNOTATION.out.annotation_ch.view { id, txt_file ->
         "SUCCESS (Annovar): Annotation for sample '$id' generated text file: $txt_file"
    }
    
	// 3. run automap for ROH detection
    DETECT_ROH(ch_input_vcf)
    DETECT_ROH.out.automap_results_ch.view{
        "SUCCESS (AutoMap): ROH detection complete. Results are in directory: automap"
    }

}

// ========================================================================================
// PROCESS DEFINITION
// ========================================================================================

process VEP_ANNOTATION {
    tag "VEP Annotation on ${id}"
    publishDir "${outdir}", mode: 'copy', pattern: "*.json"

    container 'ensemblorg/ensembl-vep:release_113.0'

    cpus 4
    memory '16 GB' // VEP might need more memory, this can be changed based on situation

    input:
    tuple val(id), path(vcf), path(outdir)

    output:
    // output JSON file
    tuple val(id), path("${id}.vep.json"), emit: vep_json_ch

    script:
	// 定義custom database的路徑，此處的clinvar是從官方ftp下載的最新版，而非cache中的
    def clinvar_path_in_cache = "${params.clinvar_custom_dir}/${params.clinvar_vcf_name}"

    """
    set -e

    vep --offline --cache --assembly GRCh37 \\
        -i ${vcf} \\
        --dir_cache /vep_cache \\
        --refseq \\
        --fasta /ref_dir/${params.ref_fasta} \\
        --no_stats \\
        --hgvs \\
        --numbers \\
        --canonical \\
        --total_length \\
        --biotype \\
        --output_file ${id}.vep.json \\
        --symbol \\
        --variant_class \\
        --json \\
        --fork ${task.cpus} \\
        --custom file=/vep_cache/${clinvar_path_in_cache},short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN%CLNSIGCONF \\
        --force_overwrite
    """
}

process DETECT_ROH{
    tag "ROH detection on ${id}"
    publishDir "${outdir}", mode: 'copy'

    container 'hydragenetics/automap:1.2'
    
    cpus 1
    memory '4 GB'

    input:
    tuple  val(id), path(vcf), path(outdir)

    output:
    // 將包含所有結果的 automap 目錄作為輸出
    path("automap"), emit: automap_results_ch

    script:
    """
    # 執行 automap 命令，並將結果輸出到一個名為 automap 的目錄中
    mkdir automap
    
    automap \\
        --vcf ${vcf} \\
        --genome hg19 \\
        --out automap \\
        --chrX
    """
}

process ANNOVAR_ANNOTATION {
    tag "Annovar annotation on ${id}"
    publishDir "$outdir", mode: 'copy', pattern: "${id}.*"

    cpus 4

    input:
    tuple val(id), path(vcf), path(outdir)

    output:
    tuple val(id), path("${id}.avinput")              , emit: avinput_ch
    tuple val(id), path("${id}.*.txt")                , emit: annotation_ch

    script:
    // 將 Annovar 的路徑定義為參數，使其更靈活
    def annovar_script_path = "${params.annovar_path}/annovar_scripts"
    def humandb_path = "${params.annovar_path}/humandb/"
    
    // 使用 id 來命名所有中間和輸出檔案，確保唯一性
    def avinput_file = "${id}.avinput"
    def annovar_prefix = "${id}" // Annovar 會自動加上 hg19_multianno.txt 等後綴

    """
    # [本地執行]：不使用容器，直接在本地執行 Perl 腳本
    set -e

    # 步驟 1: 將 VCF 轉換為 Annovar 的 avinput 格式
    perl ${annovar_script_path}/convert2annovar.pl \\
        -format vcf4 \\
        ${vcf} \\
        -withzyg \\
        -include \\
        -outfile ${avinput_file}

    # 步驟 2: 執行 table_annovar.pl 進行註解
    perl ${annovar_script_path}/table_annovar.pl \\
        ${avinput_file} \\
        ${humandb_path} \\
        -buildver hg19 \\
        --polish \\
        --intronhgvs 20 \\
        -out ${annovar_prefix} \\
        -remove \\
        -protocol refGeneWithVer,avsnp151,gff3,gnomad211_genome,twnaf_annovarin,popfreq_all_20150413,test_HGMD_annovar,LOVD_all,clinvar_20250504,intervar_20180118,dbscsnv11,spidex,dbnsfp35a \\
        -operation g,f,r,f,f,f,f,f,f,f,f,f,f \\
        --gff3dbfile hg19_rmsk.gff \\
        --argument '-hgvs,,,,,,,,,,,,' \\
        -nastring . \\
        --thread ${task.cpus} \\
        -otherinfo
    """
}
