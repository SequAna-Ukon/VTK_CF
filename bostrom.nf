#!/usr/bin/env nextflow

/*
All code associated with the analysis of the Böstrom dataset for the 2023 VTK
*/

samples_ch = Channel.fromFilePairs("/home/VTK23/sharaf/raw_reads/*/*.fastq.gz", size: 1).map{[it[0], it[1][0]]}

process fastp{
    tag "${sample}"
    conda "fastp -c bioconda"
    publishDir "/home/VTK23/Abdo/fastp/", pattern: "*.html"

    input:
    tuple val(sample), path(read_1) from samples_ch

    output:
    file "${sample}.fastp.html" into fastp_out_ch
    tuple val(sample), file("${sample}.clean.fq.gz") into kallisto_in_ch

    script:
    """
    fastp -q 20 -i $read_1 -o ${sample}.clean.fq.gz
    mv fastp.html ${sample}.fastp.html
    """
}
