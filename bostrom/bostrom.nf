#!/usr/bin/env nextflow

/*
example code associated with the analysis of the Böstrom dataset for the CF VTK
*/

nextflow.enable.dsl=2

samples_ch = Channel.fromFilePairs("/home/bio16100/bostrom/Reference/raw_reads/fastq/*.fastq.gz", size: 1).map{[it[0], it[1][0]]}

process fastp {
    tag "${sample}"
    conda "fastp -c bioconda"
    cpus 16
    publishDir "./fastp/", mode: "copy"
    input:
    tuple val(sample), path(read_1)

    output:
    file("${sample}.fastp.html")
    tuple val(sample), file("${sample}.clean.fq.gz")

    script:
    """
    fastp -q 20 --thread ${task.cpus} -i $read_1 -o ${sample}.clean.fq.gz -h ${sample}.fastp.html
    """
}

workflow {
    samples_ch.view()
    fastp(samples_ch)
}

