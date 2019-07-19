FROM ubuntu:latest
RUN rm /bin/sh && ln -s /bin/bash /bin/sh
MAINTAINER Chris Plaisier <plaisier@asu.edu> (adapted from Steve Tsang)
RUN apt-get update

RUN apt-get install --yes \
 build-essential \
 gcc-multilib \
 apt-utils \
 zlib1g-dev \
 vim-common \
 wget \
 python3 \
 python3-pip \
 git \
 pigz

# Get latest STAR source from releases
# Alternatively, get STAR source using git
RUN git clone https://github.com/alexdobin/STAR.git
WORKDIR /STAR/source/

# Build STAR
#RUN pwd
RUN make STAR

# If you have a TeX environment, you may like to build the documentation
# make manual
ENV PATH /STAR/source:$PATH
WORKDIR /

# Build index for human alignment
RUN mkdir /GRCh38.p12
WORKDIR /GRCh38.p12
RUN wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.primary_assembly.annotation.gtf.gz
RUN wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/GRCh38.primary_assembly.genome.fa.gz
RUN gunzip *.gz
RUN mkdir /index
RUN STAR --runThreadN 12 --runMode genomeGenerate --genomeDir /index --genomeFastaFiles GRCh38.primary_assembly.genome.fa --sjdbGTFfile gencode.v31.primary_assembly.annotation.gtf --sjdbOverhang 100

WORKDIR /

# Get rid of sequence information to reduce image size
RUN rm -rf /GRCh38.p12

# Install cutadapt
RUN pip3 install cutadapt

# Install ht-seq
RUN pip3 install HTSeq
