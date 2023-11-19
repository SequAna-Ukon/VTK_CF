# Use an image from the rocker/tidyverse repository
FROM gitpod/workspace-full

# Update the package list and install software using apt-get
RUN sudo apt-get update && \
    sudo apt-get install -y \
        # R environment
        r-base \
        # System libraries for R packages
        libharfbuzz-dev \
        libfribidi-dev \
        libmariadb-dev

# R Environment
RUN sudo R -e "install.packages(c('languageserver', 'tidyverse', 'BiocManager'))"
# Install Bioconductor and CRAN packages using BiocManager
RUN sudo R -e "BiocManager::install(c('RMariaDB', 'tximport', 'DESeq2', 'pheatmap', 'vsn', 'matrixStats', 'GenomicFeatures'), ask = FALSE, force = TRUE, dependencies = TRUE)"
RUN sudo R -e "install.packages('devtools')"

# Download and Install Quarto
RUN sudo curl -LO https://quarto.org/download/latest/quarto-linux-amd64.deb && \
    sudo apt-get install -y gdebi-core && \
    sudo gdebi -n quarto-linux-amd64.deb && \
    rm quarto-linux-amd64.deb
# Modify the PATH variable
ENV PATH=/opt/quarto/bin/tools:$PATH

# Copy DESeq2 R Script to /workspace
COPY kallisto_out /home/gitpod/kallisto_out

# Define the command to run when the container starts
CMD ["sudo", "-s", "/bin/bash"]
