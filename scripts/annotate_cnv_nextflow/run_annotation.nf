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
         COPY NUMBER VARIANT ANNOTATION PIPELINE
         ===================================
         Input Files      : ${params.input_vcf}
         Output Dir       : ${params.output_dir}
         ---
         Execution reports will be generated in '${params.output_dir}'.
         Check nextflow.config for trace/report settings.
         """
		 

// ========================================================================================
// WORKFLOW DEFINITIONS
// ========================================================================================

workflow {
    Channel
        .fromPath(params.input_vcf)
        .map { vcf_file ->
            def id = vcf_file.simpleName
            return [id, vcf_file] 
        }
        .view()
        .set { ch_input_vcf }

    ch_annotsv_result = ANNOTSV_ANNOTATION( ch_input_vcf)
    ch_annotsv_result.view{ id, tsv ->
        "success: annotate for sample '$id' generated $tsv"
    }
    GENERATE_HTML_REPORT( ch_annotsv_result)
    GENERATE_HTML_REPORT.out.knot_html_ch.view {id, html_file->
        "SUCCESS (AnnotSV): Annotation for sample '$id' generated html file: $html_file"  
    }
}

// =========================================================================================
// PROCESS DEFINITIONS
// =========================================================================================

process ANNOTSV_ANNOTATION {
    tag "AnnotSV Annotation on ${id}"
    publishDir "${params.output_dir}", mode: "copy"

    container 'quay.io/biocontainers/annotsv:3.4.6--py313hdfd78af_0'
  
    cpus 4
    memory '4 GB'
  
    input:
    tuple val(id), path(vcf_file)
  
    output:
    tuple val(id), path("annotSV/${id}_annotsv.tsv"), emit: annotsv_results_ch
  
    script:
    def annotation_dir="/annotsv_bundle/AnnotSV_annotations/"
    """
    set -e
    
    mkdir annotSV

    AnnotSV -SVinputFile ${vcf_file} \\
            -genomeBuild GRCh37 \\
            -outputDir annotSV \\
            -outputFile annotSV/${id}_annotsv.tsv \\
            -annotationsDir ${annotation_dir}
    """
}


process GENERATE_HTML_REPORT {
    tag "Generate knotAnnotSV html report on ${id}"
    publishDir "${params.output_dir}", mode: "copy"
	
    container "knotannotsv"
    cpus 1
	
    input:
    tuple val(id), path(annotsv_tsv)
	
    output:
    tuple val(id), path("annotSV/knot_${id}_annotsv.html"), emit: knot_html_ch
	
    script:
    def yaml="/annotsv_bundle/config_AnnotSV.yaml"
    """
    set -e

    mkdir annotSV

    knotAnnotSV --annotSVfile ${annotsv_tsv} \\
        --configFile ${yaml} \\
        --outDir annotSV \\
        --outPrefix knot
    """

}


