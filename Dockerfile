FROM ubuntu:22.04

################################################################################
# Genotype Imputation Pipeline Docker Image
#
# This Docker image contains all the tools required to run the imputation
# pipeline: PLINK, PLINK2, bcftools, vcftools, samtools, tabix, R, Java,
# Python, GenotypeHarmonizer, ADMIXTURE, Eagle, and Minimac4.
################################################################################

# Avoid interactive prompts during installation
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC

# Install system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    git \
    unzip \
    bzip2 \
    build-essential \
    autoconf \
    automake \
    libtool \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libncurses5-dev \
    parallel \
    python3 \
    python3-pip \
    openjdk-11-jdk \
    r-base \
    r-base-dev \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /opt

################################################################################
# Install bcftools, samtools, and tabix (htslib)
################################################################################

RUN wget https://github.com/samtools/htslib/releases/download/1.18/htslib-1.18.tar.bz2 && \
    tar -xjf htslib-1.18.tar.bz2 && \
    cd htslib-1.18 && \
    ./configure --prefix=/usr/local && \
    make && make install && \
    cd .. && rm -rf htslib-1.18 htslib-1.18.tar.bz2

RUN wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2 && \
    tar -xjf samtools-1.18.tar.bz2 && \
    cd samtools-1.18 && \
    ./configure --prefix=/usr/local && \
    make && make install && \
    cd .. && rm -rf samtools-1.18 samtools-1.18.tar.bz2

RUN wget https://github.com/samtools/bcftools/releases/download/1.18/bcftools-1.18.tar.bz2 && \
    tar -xjf bcftools-1.18.tar.bz2 && \
    cd bcftools-1.18 && \
    ./configure --prefix=/usr/local && \
    make && make install && \
    cd .. && rm -rf bcftools-1.18 bcftools-1.18.tar.bz2

################################################################################
# Install vcftools
################################################################################

RUN wget https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz && \
    tar -xzf vcftools-0.1.16.tar.gz && \
    cd vcftools-0.1.16 && \
    ./configure --prefix=/usr/local && \
    make && make install && \
    cd .. && rm -rf vcftools-0.1.16 vcftools-0.1.16.tar.gz

################################################################################
# Copy pre-existing tools from required_tools/
################################################################################

# Copy PLINK 1.9
COPY required_tools/plink /usr/local/bin/plink
RUN chmod +x /usr/local/bin/plink

# Copy PLINK 2.0
COPY required_tools/plink2 /usr/local/bin/plink2
RUN chmod +x /usr/local/bin/plink2

# Copy Eagle v2.4.1
COPY required_tools/Eagle_v2.4.1 /opt/eagle
RUN ln -s /opt/eagle/eagle /usr/local/bin/eagle && \
    chmod +x /opt/eagle/eagle

# Copy Minimac4 v4.1.6
COPY required_tools/minimac4 /usr/local/bin/minimac4
RUN chmod +x /usr/local/bin/minimac4

################################################################################
# Install ADMIXTURE
################################################################################

RUN wget http://dalexander.github.io/admixture/binaries/admixture_linux-1.3.0.tar.gz && \
    tar -xzf admixture_linux-1.3.0.tar.gz && \
    mv dist/admixture_linux-1.3.0/admixture /usr/local/bin/ && \
    chmod +x /usr/local/bin/admixture && \
    rm -rf dist admixture_linux-1.3.0.tar.gz

################################################################################
# Copy GenotypeHarmonizer
################################################################################

COPY required_tools/GenotypeHarmonizer /opt/GenotypeHarmonizer

################################################################################
# Install liftOver tool
################################################################################

RUN wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver && \
    chmod +x liftOver && \
    mv liftOver /usr/local/bin/

################################################################################
# Install Python packages
################################################################################

RUN pip3 install --no-cache-dir numpy pandas scipy matplotlib seaborn

################################################################################
# Install R packages
################################################################################

RUN R -e "install.packages(c('data.table', 'ggplot2', 'dplyr', 'tidyr'), repos='https://cloud.r-project.org/')"

################################################################################
# Copy required_tools directory and pipeline scripts
################################################################################

COPY required_tools /pipeline/required_tools
COPY scripts/ /pipeline/scripts/
COPY imputation_pipeline.sh /pipeline/imputation_pipeline.sh
COPY prepare_reference_panel.sh /pipeline/prepare_reference_panel.sh

RUN chmod +x /pipeline/imputation_pipeline.sh && \
    chmod +x /pipeline/prepare_reference_panel.sh

################################################################################
# Set up environment
################################################################################

ENV PATH="/usr/local/bin:${PATH}"
ENV JAVA_HOME="/usr/lib/jvm/java-11-openjdk-amd64"

# Set working directory for pipeline
WORKDIR /pipeline

################################################################################
# Re-install Minimac4 to ensure correct setup
################################################################################
# RUN wget https://github.com/statgen/Minimac4/releases/download/v4.1.6/minimac4-4.1.6-Linux-x86_64.sh
# RUN yes | bash minimac4-4.1.6-Linux-x86_64.sh -b -p /usr/local/bin/minimac4 > /mynimac4_install.log 2>&1
# RUN mv minimac4-4.1.6-Linux-x86_64 minimac4
# RUN mv /pipeline/minimac4/bin/minimac4 /usr/local/bin/minimac4
# RUN chmod +x /usr/local/bin/minimac4
# RUN rm minimac4-4.1.6-Linux-x86_64.sh
# RUN rm -rf /pipeline/minimac4

################################################################################
# Verify installations
################################################################################

RUN echo "Verifying installations..."
RUN bcftools --version
RUN samtools --version
RUN vcftools --version
RUN plink --version
RUN plink2 --version
RUN minimac4 --version
RUN admixture --version
RUN which liftOver
RUN java -version
RUN python3 --version
RUN R --version
RUN parallel --version
RUN echo "All tools installed successfully!"

################################################################################
# Entry point
################################################################################

CMD ["/bin/bash"]
