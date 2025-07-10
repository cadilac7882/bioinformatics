---

# 🧬 Secondary Analysis Pipeline (GPU-based, Dockerized)

This repository includes scripts for performing secondary analysis (from FASTQ to VCF), quality control, and downstream genomic evaluations. The pipeline is designed to run in a Docker environment and is GPU-accelerated (Parabricks), where applicable.

---

## 📁 Pipeline Overview

### 🔄 `automation_full.sh`

* **Purpose**:
  Full automation of the secondary analysis pipeline. Includes:

  * Alignment
  * Variant calling
  * Basic quality control metrics
* **Notes**: GPU-based; requires Parabricks

---

## 📁 FASTQ Management

### 📦 `archive_fastq.sh`

* **Purpose**:
  Archive `*.fastq.gz` files located in the `fastq/` folder based on `sampleList.txt`.

---

## 🧪 Mapping Quality Evaluation

### 🧾 `auto_mapq_docker_ver.sh` *(wrapper)*

### 🧾 `testMapQuality.sh` *(single-thread)*

### 🧾 `testMapQuality_multithreads.sh` *(multi-threaded)*

* **Purpose**:
  Estimate mapping quality in **target regions** for each BAM file.

---

## 📊 BAM Quality Metrics

### 🧾 `auto_QC_metric_docker_ver.sh` *(wrapper)*

### 🧾 `samtools_metric_multithreads.sh` *(multi-threaded)*

* **Purpose**:
  Compute alignment statistics (e.g., coverage, duplication rate) using `samtools` and `bed` files.

---

## ⚙️ Secondary Analysis with Parabricks

### 🧾 `nv_run_docker_ver.sh`

* **Purpose**: Germline variant calling on **hg19**

### 🧾 `nv_run_docker_ver_hg38.sh`

* **Purpose**: Germline variant calling on **hg38**

### 🧾 `nv_run_docker_ver_somatic.sh`

* **Purpose**: Somatic variant calling pipeline
* **Toolset**: GPU-accelerated Parabricks (bwa-mem2, HaplotypeCaller/Mutect2)

---

## 🧼 Variant Filtration

### 🧾 `hardfiltering_docker_ver.sh`

* **Purpose**:
  Apply **hard-filtering** to VCF files using GATK best-practices.

---

## 🧬 Copy Number Variation (CNV)

### 🧾 `test_gatk_germlineCNV.sh`

* **Purpose**:
  Run GATK Germline CNV caller in **cohort mode**.

### 🧾 `merge_contig_ploidy.sh`

* **Purpose**:
  Merge contig-level ploidy results predicted by GATK gCNV model.

---

## 🧬 Exon-Level Coverage

### 🧾 `mosdetph_cds.sh`

* **Purpose**:
  Compute exon-level depth metrics using `mosdepth` and a BED file of coding exons.

---

## ✂️ Preprocessing

### 🧾 `trim.sh`

* **Purpose**:
  Trim raw FASTQ files using **Trimmomatic**.

---

## 🔗 Phasing

### 🧾 `whatshap.sh`

* **Purpose**:
  Perform **read-based phasing** using WhatsHap.

---

## 📌 Notes

* All scripts are designed to work with Docker and/or Parabricks environments.
* Ensure proper volume mapping and input paths when executing scripts.
* Sample names should match across FASTQ, BAM, and sample list files.

---
