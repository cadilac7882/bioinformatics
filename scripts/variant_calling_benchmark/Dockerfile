FROM ubuntu:18.04

ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies (based on official dockerfile)
RUN apt-get update && \
    apt-get install -y \
        autoconf \
        build-essential \
        bzip2 \
        cmake \
        cython \
        git \
        libbz2-dev \
        libncurses5-dev \
        openjdk-8-jdk \
        pkg-config \
        python \
        python2.7 \
        python2.7-dev \
        python-setuptools \
        python-pip \
        python-psutil \
        python-numpy \
        python-pandas \
        python-distribute \
        python-pysam \
        python-scipy \
        software-properties-common \
        wget \
        zlib1g-dev && \
    apt-get clean -y

# Fix the bx-python Unicode issue by building from source
RUN pip install --no-binary bx-python bx-python

# Clone hap.py repository
RUN git clone https://github.com/Illumina/hap.py.git /opt/hap.py-source

# Create data directory
RUN mkdir -p /opt/hap.py-data

# Install Apache Ant (required for RTG tools)
WORKDIR /opt
RUN wget http://archive.apache.org/dist/ant/binaries/apache-ant-1.9.7-bin.tar.gz && \
    tar xzf apache-ant-1.9.7-bin.tar.gz && \
    rm apache-ant-1.9.7-bin.tar.gz

ENV PATH $PATH:/opt/apache-ant-1.9.7/bin

# Install hap.py
WORKDIR /opt/hap.py-source
RUN python install.py /opt/hap.py --with-rtgtools --no-tests

# Set working directory
WORKDIR /opt/hap.py

# Run basic tests
RUN bin/test_haplotypes

# Clean up source folder
RUN rm -rf /opt/hap.py-source

# Set up environment
ENV PATH="/opt/hap.py/bin:${PATH}"

# Set working directory for data
WORKDIR /data

# Default command
CMD ["hap.py", "--help"]
