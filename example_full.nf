nextflow.enable.dsl=2

samples_ch = Channel.fromFilePairs("./raw_reads/fastq/*.fastq.gz", size: 1).map { [it[0], it[1][0]] }

process fastp {
    tag "${sample}"
    conda "fastp -c bioconda"
    cpus 16
    publishDir "./fastp/", mode: "copy"

    input:
    tuple val(sample), path(read_1)

    output:
    file("${sample}.fastp.html")
    file("${sample}.clean.fq.gz")


    script:
    """
    fastp -q 20 --thread ${task.cpus} -i $read_1 -o ${sample}.clean.fq.gz -h fastp.html ${sample}.fastp.html
    """
}

process kallistoIndex {
    tag "kallistoIndex"
    conda 'kallisto -c bioconda'

    publishDir path: './results/', mode: 'copy'

    input:
    path transcriptome_fa

    output:
    file "kallisto_index"

    script:
    """
    kallisto index -i kallisto_index ${transcriptome_fa}
    """
}


workflow {
    input_ch = samples_ch
    transcriptome_fa = file("./Reference/transcriptome/Homo_sapiens.GRCh38.cdna.all.fa.gz")

    // Invoke fastp before kallisto
    fastp(input_ch)

    // Invoke kallistoIndex before kallisto
    kallistoIndex(transcriptome_fa)
}
