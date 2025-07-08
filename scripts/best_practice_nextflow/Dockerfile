# Dockerfile to install BWA, Samtools, and GATK using Mamba (Multi-stage build)

# --- STAGE 1: The Builder ---
# 使用 Mambaforge 作為基礎來安裝所有軟體
FROM condaforge/mambaforge:23.11.0-0 AS builder

# 從 bioconda 和 conda-forge channel 安裝所有工具的指定版本
# Mamba 會自動解析並安裝所有依賴，包括 GATK 所需的 Java 8
RUN mamba install -y -c bioconda -c conda-forge \
        bwa=0.7.17 \
        samtools=1.19 \
        gatk4=4.3.0.0 \
        bcftools=1.21 \
        bedtools=2.31.0 && \
    # 清理快取以減小最終複製的檔案大小
    mamba clean -afy

# --- STAGE 2: The Final Image ---
# 從一個非常小、乾淨的基礎映像開始
FROM debian:bullseye-slim

# ==============================================================================
# **[THE FIX]**: Install the 'procps' package.
# This provides the `ps` command, which Nextflow requires for task monitoring
# and generating execution trace files. We use --no-install-recommends to
# keep the installation minimal.
RUN apt-get update && \
    apt-get install -y --no-install-recommends procps && \
    # Clean up apt cache to keep the final image small.
    rm -rf /var/lib/apt/lists/*
# ==============================================================================

# 從 builder 階段僅複製已安裝好的 Conda 環境到最終映像中
COPY --from=builder /opt/conda /opt/conda

# 將 Conda 的 bin 目錄添加到系統的 PATH 環境變數中
# 這樣我們就可以直接呼叫 bwa, samtools, gatk
ENV PATH=/opt/conda/bin:$PATH

# 設置工作目錄
WORKDIR /workdir

# **[重要]**：不設置 ENTRYPOINT 或 CMD，讓 Nextflow 完全控制命令的執行
