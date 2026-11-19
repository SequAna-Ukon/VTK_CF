/*
example code associated with the analysis of the Böstrom dataset for the CF VTK
*/

nextflow.enable.dsl=2

// Input FASTQ files
params.fastq_files = "/home/bio16100/bostrom/Reference/raw_reads/fastq/*.fastq.gz"

process fastp {
    tag "${sample}"
    conda "fastp -c bioconda"
    cpus 10
    publishDir "./fastp/", mode: "copy"

    input:
    tuple val(sample), path(fastq)

    output:
    file("${sample}.fastp.html")
    tuple val(sample), file("${sample}.clean.fq.gz")

    script:
    """
    fastp -q 20 --thread ${task.cpus} -i $fastq -o ${sample}.clean.fq.gz -h ${sample}.fastp.html
    """
}

workflow {
    samples_ch = Channel.fromFilePairs(params.fastq_files, size: 1).map{[it[0], it[1][0]]}

    fastp(samples_ch)
}
