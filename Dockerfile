# Use an image from the rocker/tidyverse repository
FROM rocker/tidyverse:latest

# Update the package list and install software using apt-get
RUN apt-get update && \
    apt-get install -y \
        less \
        libghc-bzlib-dev \
        curl \
        python3-pip \
        openjdk-11-jdk \
        fastp \
        kallisto

# Install the 'languageserver' package using R
RUN R -e "install.packages('languageserver')"

# Install the 'BiocManager' package using R
RUN R -e "install.packages('BiocManager')"

# Install Bioconductor packages using BiocManager
RUN R -e "BiocManager::install(c('RMariaDB', 'tximport', 'DESeq2', 'pheatmap', 'vsn', 'matrixStats', 'GenomicFeatures'), ask = FALSE, force = TRUE, dependencies = TRUE)"

# Define the command to run when the container starts
CMD ["/bin/bash"]
