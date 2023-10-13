FROM rocker/tidyverse:latest

RUN apt-get update && \
    apt-get install -y \
        fastp \
        kallisto

# For Quarto VS Code Extention
RUN R -e "install.packages('languageserver')"

# For Bioconductor
RUN R -e "install.packages('BiocManager')"

# For DiffExpr Analysis
RUN R -e "BiocManager::install(c('RMariaDB', 'tximport', 'DESeq2', 'pheatmap', 'vsn', 'matrixStats', 'GenomicFeatures'))"

CMD ["/bin/bash"]