# Use an official Ubuntu base image
FROM ubuntu:20.04

# Maintainer information
MAINTAINER Manish Kumar <mandecent.gupta@gmail.com>

# Define arguments for versions
ARG star_version
ARG kallisto_version
ARG samtools_version
ARG container_timezone=Etc/UTC

# Set timezone
RUN ln -snf /usr/share/zoneinfo/$container_timezone /etc/localtime && \
    echo $container_timezone > /etc/timezone

# Install basic dependencies and tools
RUN apt-get update && \
    apt-get install -y wget skewer git build-essential \
                       libncurses5-dev zlib1g-dev libbz2-dev \
                       liblzma-dev libcurl3-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Create directory for software installations
RUN mkdir -p /software/utils

# Download and install Kallisto
WORKDIR /tmp
RUN wget https://github.com/pachterlab/kallisto/releases/download/$kallisto_version/kallisto_linux-$kallisto_version.tar.gz && \
    tar -xzf kallisto_linux-$kallisto_version.tar.gz -C /tmp && \
    mv /tmp/kallisto /software/utils

# Download and install STAR
RUN wget https://github.com/alexdobin/STAR/archive/$star_version.tar.gz && \
    tar -xzf $star_version.tar.gz -C /tmp && \
    mv /tmp/STAR-$star_version/bin/Linux_x86_64_static/STAR /software/utils

# Download, compile, and install Samtools
WORKDIR /tmp
RUN wget https://github.com/samtools/samtools/releases/download/$samtools_version/samtools-$samtools_version.tar.bz2 && \
    tar jxf samtools-$samtools_version.tar.bz2 && \
    rm samtools-$samtools_version.tar.bz2 && \
    cd samtools-$samtools_version && \
    ./configure --prefix=/software/utils && \
    make && \
    make install

# Clean up temporary files
RUN rm -rf /tmp/*

# Set the working directory
WORKDIR /software/utils

# Set entrypoint (optional)
# ENTRYPOINT ["/bin/bash"]