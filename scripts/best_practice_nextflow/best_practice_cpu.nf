#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """
         G A T K   B E S T   P R A C T I C E S   W O R K F L O W
         ===================================
         Input Files      : ${params.input}
         Reference Dir    : ${params.ref_dir} (mounted as /refdir in container)
         Output Dir       : ${params.outdir}
         ---
         Execution reports will be generated in '${params.outdir}'.
         Check nextflow.config for trace/report settings.
         """

workflow {
    // 1. Create a channel for paired-end FASTQ files
    Channel
        .fromFilePairs(params.input)
        .ifEmpty { exit 1, "Cannot find any FASTQ files matching: ${params.input}" }
        .set { ch_fastq_pairs }

    // 2. Pre-processing, BQSR, and Variant Calling chain
    ch_sorted_bams = ALIGN_AND_SORT(ch_fastq_pairs)
    MARK_DUP(ch_sorted_bams)

    ch_marked_dup_bam = MARK_DUP.out.marked_bam_ch
    ch_with_bqsr_report = BQSR(ch_marked_dup_bam)
    ch_final_bam = APPLY_BQSR(ch_with_bqsr_report)

    ch_raw_vcf = HAPLOTYPE_CALLER(ch_final_bam)

    // 3. Post-processing and hard-filtering chain
    ch_postproc_vcf = POSTPROC_VCF(ch_raw_vcf)
    HARDFILTERING(ch_postproc_vcf)

    // 4. View final output for debugging
    HARDFILTERING.out.hf_vcf_ch.view { id, vcf ->
        "Final hard-filtered VCF for sample '$id' is at: $vcf"
    }
}

// ========================================================================================
// PROCESS DEFINITIONS
// ========================================================================================

process ALIGN_AND_SORT {
    tag "Align and sort BAM on ${id}"

    container 'dgm/secondary_analysis_toolkits:1.0'
    cpus 32
    memory '32 GB'
    
    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("${id}.sorted.bam")

    script:
    def (r1, r2) = reads
    def ref_path_in_container = "/refdir/${params.ref_fasta}"
    """
    bwa mem \\
        -t ${task.cpus} \\
        -K ${params.bwa_k_option} \\
        -R '@RG\\tID:${id}\\tLB:lib1\\tPL:bar\\tSM:${id}\\tPU:${id}' \\
        ${ref_path_in_container} \\
        ${r1} \\
        ${r2} | \\
    gatk SortSam \\
        --java-options "-Xmx30g" \\
        -I /dev/stdin \\
        -O ${id}.sorted.bam \\
        --SORT_ORDER coordinate
    """
}

process MARK_DUP {
    tag "Mark Duplicates on ${id}"
    // Publish to the same central directory, using 'pattern' to be specific
    publishDir "${params.outdir}/${id}", mode: 'copy', pattern: "*.metrics.txt"

    container 'dgm/secondary_analysis_toolkits:1.0'

    // MarkDuplicates is more memory/IO bound than CPU bound.
    // Requesting fewer CPUs but ample memory is more efficient.
    cpus 8
    memory '32 GB'

    input:
    tuple val(id), path(bam)

    output:
    // Capture BOTH output files.
    // The 'emit' keyword lets us name the output channels for clarity.
    tuple val(id), path("${id}.sorted.markdup.bam"), emit: marked_bam_ch
    path("${id}.metrics.txt")                     , emit: metrics_ch

    script:
    // Name the metrics file with the sample ID to keep it unique.
    """
    gatk MarkDuplicates \\
      --java-options "-Xmx30g" \\
      -I ${bam} \\
      -O ${id}.sorted.markdup.bam \\
      -M ${id}.metrics.txt
    """
}

process BQSR {
    tag "Generate BQSR Report on ${id}"
    publishDir "${params.outdir}/${id}", mode: 'copy', pattern: "*.txt"

    container 'dgm/secondary_analysis_toolkits:1.0'
    cpus 8
    memory '32 GB'
    
    input:
    tuple val(id), path(bam)

    output:
    tuple val(id), path(bam), path("${id}.bqsr_recal.txt"), emit: bqsr_report_ch

    script:
    def ref_path_in_container = "/refdir/${params.ref_fasta}"
    def dbsnp_sites = "/refdir/dbsnp_138.hg19.vcf.gz"
    def indel_sites_1000G = "/refdir/1000G_phase1.indels.hg19.sites.vcf.gz"
    def indel_sites_mills = "/refdir/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"
	
    """
    gatk BaseRecalibrator \\
      --java-options "-Xmx30g" \\
      --input ${bam} \\
      --output ${id}.bqsr_recal.txt \\
      --known-sites ${dbsnp_sites} \\
      --known-sites ${indel_sites_1000G} \\
      --known-sites ${indel_sites_mills} \\
      --reference ${ref_path_in_container}
    """
}

process APPLY_BQSR {
    tag "Apply BQSR on ${id}"
    publishDir "${params.outdir}/${id}", mode: 'copy', pattern: "*sorted.markdup.bqsr.bam*"

    container 'dgm/secondary_analysis_toolkits:1.0'
    cpus 8
    memory '32 GB'
    
    input:
    tuple val(id), path(bam), path(bqsr_report)

    output:
    tuple val(id), path("${id}.sorted.markdup.bqsr.bam"), emit: bqsr_bam_ch

    script:
    def ref_path_in_container = "/refdir/${params.ref_fasta}"
	
    """
    gatk ApplyBQSR \\
      --java-options "-Xmx30g" \
      -R ${ref_path_in_container} \\
      -I ${bam} \\
      --bqsr-recal-file ${bqsr_report} \\
      -O ${id}.sorted.markdup.bqsr.bam
    """
}

process HAPLOTYPE_CALLER {
    tag "Haplotype Caller on ${id}"
    publishDir "${params.outdir}/${id}", mode: 'copy', pattern: "*.hc.vcf*"

    container 'dgm/secondary_analysis_toolkits:1.0'
    cpus 32
    memory '32 GB'
	
    
    input:
    tuple val(id), path(bam)

    output:
    tuple val(id), path("${id}.hc.vcf"), emit: hc_vcf_ch

    script:
    def ref_path_in_container = "/refdir/${params.ref_fasta}"
	
    """
    gatk HaplotypeCaller \\
      --java-options "-Xmx30g" \\
      --input ${bam} \\
      --output ${id}.hc.vcf \\
      --reference ${ref_path_in_container} \\
      --native-pair-hmm-threads ${task.cpus}
    """
}

process POSTPROC_VCF {
    tag "POST-processing VCF on ${id}"
	
    container 'dgm/secondary_analysis_toolkits:1.0'
    cpus 4
	
    input:
    tuple val(id), path(vcf)
	
    output:
    tuple val(id), path("${id}.hc.postproc.vcf"), emit: postproc_vcf_ch
	
    script:
    def ref_path_in_container = "/refdir/${params.ref_fasta}"
	
    """
    bcftools +fill-tags ${vcf} -- -t VAF | \\
    bcftools norm -m-any - | \\
    bcftools norm -f ${ref_path_in_container} - -o ${id}.hc.postproc.vcf
    """	
}

process HARDFILTERING {
    tag "Hard-filtering VCF on ${id}"
    publishDir "${params.outdir}/${id}", mode: 'copy', pattern: "*.hc.postproc.hf.vcf*"
	
    container 'dgm/secondary_analysis_toolkits:1.0'
    cpus 4
    memory '8 GB'
	
    input:
    tuple val(id), path(vcf)
	
    output:
    tuple val(id), path("${id}.hc.postproc.hf.vcf"), emit: hf_vcf_ch
	
    script:
	
    """
    # [優化] 使用 set -e 使腳本更簡潔、更健壯
    set -e
	
    # 1. 將 SNP 和 INDEL 分離到不同的 VCF 檔案
    gatk SelectVariants \\
      -V ${vcf} \\
      -O ${id}.snp.vcf \\
      -select-type SNP
	  
    gatk SelectVariants \\
      -V ${vcf} \\
      -O ${id}.indel.vcf \\
      -select-type INDEL
	
    # 2. 分別對 SNP 和 INDEL 應用硬過濾(hard filtering)
    gatk VariantFiltration \\
      -V ${id}.snp.vcf \\
      -O ${id}.snp.filtered.vcf \\
      -filter "QD < 2.0" --filter-name "QD2" \\
      -filter "QUAL < 30.0" --filter-name "QUAL30" \\
      -filter "SOR > 3.0" --filter-name "SOR3" \\
      -filter "FS > 60.0" --filter-name "FS60" \\
      -filter "MQ < 40.0" --filter-name "MQ40" \\
      -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \\
      -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"

    gatk VariantFiltration \\
      -V ${id}.indel.vcf \\
      -O ${id}.indel.filtered.vcf \\
      -filter "QD < 2.0" --filter-name "QD2" \\
      -filter "QUAL < 30.0" --filter-name "QUAL30" \\
      -filter "FS > 200.0" --filter-name "FS200" \\
      -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"
      
    # 3. 將過濾後的 SNP 和 INDEL 合併回一個最終的 VCF 檔案
    gatk MergeVcfs \\
      -I ${id}.snp.filtered.vcf \\
      -I ${id}.indel.filtered.vcf \\
      -O ${id}.hc.postproc.hf.vcf
    """
}
