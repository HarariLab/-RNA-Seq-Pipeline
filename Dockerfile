##########################################################################################
# Dockerfile for combined RNAseq pipelinev1; modified from gtex-pipeline/rnaseq/Dockerfile
##########################################################################################

# Pull the base image
FROM ubuntu:18.04

MAINTAINER Abdallah Eteleeb <eteleeb@wustl.edu>

LABEL Image for RNAseq pre-processing and alignment pipeline

ENV DEBIAN_FRONTEND=nonintercative

RUN apt-get update && apt-get install -y software-properties-common && add-apt-repository -y ppa:openjdk-r/ppa && \
    apt-get update && apt-get install -y \
        autoconf \
        build-essential \
        cmake \
        curl \
        default-jre \
        git \
        gcc \
        g++ \
        gfortran \
        libnss-sss \
        libboost-all-dev \
        lbzip2 \
        libbz2-dev \
        libcurl3-dev \
	libhdf5-dev \
        liblzma-dev \
        libncurses5-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        libudunits2-dev \
        liblzo2-dev \
        make \
        #openjdk-7-jdk \
        openjdk-8-jdk \
        perl \
        pbzip2 \
        pigz \
        aria2 \
        python3 \
        python3-pip \
        rsync \
        unzip \
        vim-common \
        wget \
        ca-certificates \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*


#-----------------------------
# Pipeline components
#-----------------------------

# htslib (required for samtools) - Updated to 1.9
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && \
    tar -xf htslib-1.9.tar.bz2 && rm htslib-1.9.tar.bz2 && cd htslib-1.9 && \
    ./configure --enable-libcurl --enable-s3 --enable-plugins --enable-gcs && \
    make && make install && make clean


# samtools - Updated to 1.9
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    tar -xf samtools-1.9.tar.bz2 && rm samtools-1.9.tar.bz2 && cd samtools-1.9 && \
    ./configure --with-htslib=/opt/htslib-1.9 && make && make install && make clean

# bamtools
RUN cd /opt && \
    wget --no-check-certificate https://github.com/pezmaster31/bamtools/archive/v2.4.1.tar.gz && \
    tar -xf v2.4.1.tar.gz && rm v2.4.1.tar.gz && cd bamtools-2.4.1 && mkdir build && cd build && cmake .. && make && make install && make clean
ENV LD_LIBRARY_PATH /usr/local/lib/bamtools:$LD_LIBRARY_PATH

# Picard tools - Updated to 2.18.29 (so can handle SE reads)
RUN mkdir /opt/picard-tools && \
    wget --no-check-certificate -P /opt/picard-tools/ https://github.com/broadinstitute/picard/releases/download/2.18.29/picard.jar

# STAR v2.5.3a - Updated to 2.7.1a
RUN cd /opt && \
    wget --no-check-certificate https://github.com/alexdobin/STAR/archive/2.7.1a.tar.gz && \
    tar -xf 2.7.1a.tar.gz && rm 2.7.1a.tar.gz && \
    make STAR -C STAR-2.7.1a/source && make STARlong -C STAR-2.7.1a/source && \
    mv STAR-2.7.1a/source/STAR* STAR-2.7.1a/bin/Linux_x86_64/
ENV PATH /opt/STAR-2.7.1a/bin/Linux_x86_64:$PATH

# RSEM v1.3.0
RUN cd /opt && \
    wget --no-check-certificate https://github.com/deweylab/RSEM/archive/v1.3.0.tar.gz && \
    tar -xvf v1.3.0.tar.gz && rm v1.3.0.tar.gz && cd RSEM-1.3.0 && make
ENV PATH /opt/RSEM-1.3.0:$PATH

# RNA-SeQC
RUN cd /opt && \
    wget --no-check-certificate https://github.com/francois-a/rnaseqc/releases/download/v1.1.9/RNA-SeQC_1.1.9.zip && \
    unzip RNA-SeQC_1.1.9.zip -d RNA-SeQC_1.1.9 && rm RNA-SeQC_1.1.9.zip

# python modules
RUN pip3 install --upgrade pip
RUN pip3 install tables numpy pandas feather-format simplejson pytz
# numpy dependencies:
RUN pip3 install pyBigWig

# Install RSeQC
RUN pip install RSeQC

# kallisto v0.43.1
RUN cd /opt && \
    wget https://github.com/pachterlab/kallisto/releases/download/v0.43.1/kallisto_linux-v0.43.1.tar.gz && \
    tar -xf kallisto_linux-v0.43.1.tar.gz && rm kallisto_linux-v0.43.1.tar.gz
ENV PATH $PATH:/opt/kallisto_linux-v0.43.1

# bedtools
RUN cd /opt && \
    wget --no-check-certificate https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz && \
    tar -xf bedtools-2.27.1.tar.gz && rm bedtools-2.27.1.tar.gz && \
    cd bedtools2 && make && make install && make clean

# UCSC tools
RUN mkdir /opt/ucsc && \
    wget --no-check-certificate -P /opt/ucsc/ http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph && \
    wget --no-check-certificate -P /opt/ucsc/ http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig && \
    chmod 755 /opt/ucsc/*
ENV PATH /opt/ucsc:$PATH

# gcloud
RUN export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
    echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    apt-get update -y && apt-get install google-cloud-sdk -y

# fastqc v0.11.7
RUN mkdir -p /opt/tools && \
    wget -c http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && \
    unzip fastqc_v0.11.9.zip && \
    cd FastQC && \
    chmod +x fastqc && \
    cp fastqc /usr/local/bin
ENV PATH /opt/tools:$PATH

# Install MultiQC
RUN pip install git+git://github.com/ewels/MultiQC.git
#RUN pip install git+https://github.com/ewels/MultiQC.git
#RUN pip install multiqc

# salmon v1.2.0 (updated)
RUN cd /opt && \
	wget https://github.com/COMBINE-lab/salmon/releases/download/v1.2.0/salmon-1.2.0_linux_x86_64.tar.gz && \
	tar -xf salmon-1.2.0_linux_x86_64.tar.gz && rm salmon-1.2.0_linux_x86_64.tar.gz

ENV PATH $PATH:/opt/salmon-1.2.0_linux_x86_64

# Install R 
RUN apt-get install -y r-base

# install R Packages
RUN R -e "install.packages(c('stringr'), repos='http://cran.us.r-project.org')"

# install rclone 
RUN curl https://rclone.org/install.sh | bash

# install synapseclient 
RUN pip install --upgrade synapseclient

# clean up
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get autoclean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/

#COPY hg19_GencodeCompV19.bed /usr/bin/hg19_GencodeCompV19.bed
COPY ./src /src

RUN chmod +rx /src/*.py && \
    chmod +rx /src/*.pl
ENV PATH /src:$PATH

#WORKDIR /usr/bin
RUN cp /src/pipeline_unified.py /usr/bin/ngi-rnaseq

CMD ["/bin/bash"]
