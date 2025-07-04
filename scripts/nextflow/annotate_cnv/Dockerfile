FROM ubuntu:22.04

LABEL maintainer="you@example.com"
ENV DEBIAN_FRONTEND=noninteractive

# 安裝系統與 Perl 基本工具
RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    git \
    perl \
    cpanminus \
    unzip \
    wget \
    libdbi-perl \
    libdbd-sqlite3-perl \
    libjson-perl \
    libxml-simple-perl \
    libtext-csv-perl \
    libfile-slurp-perl

# 補充常缺的 Perl 模組（來自 CPAN）
RUN cpanm YAML::XS Sort::Key::Natural

# 安裝 AnnotSV
RUN wget https://github.com/lgmgeo/AnnotSV/archive/refs/tags/v3.4.6.zip && \
    unzip v3.4.6.zip && \
    mv AnnotSV-3.4.6 /opt/AnnotSV && \
    ln -s /opt/AnnotSV/AnnotSV /usr/local/bin/AnnotSV

# 安裝 knotAnnotSV
RUN git clone https://github.com/mobidic/knotAnnotSV.git /opt/knotAnnotSV && \
    chmod +x /opt/knotAnnotSV/knotAnnotSV.pl && \
    ln -s /opt/knotAnnotSV/knotAnnotSV.pl /usr/local/bin/knotAnnotSV

WORKDIR /data
CMD ["perl", "/opt/knotAnnotSV/knotAnnotSV.pl"]
