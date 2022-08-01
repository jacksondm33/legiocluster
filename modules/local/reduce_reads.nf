process REDUCE_READS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::multiqc=1.12' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.12--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)
    val random
    val k

    output:
    tuple val(meta), path("${prefix}_reduced_[12].fastq"), emit: reduced_reads
    tuple val(meta), path(log_file)                      , emit: log
    path  "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    log_level     = "INFO"
    reduced_reads = "${prefix}_reduced_1.fastq ${prefix}_reduced_2.fastq"
    log_file      = "${prefix}.log"

    template 'reduce_reads.py'
}
