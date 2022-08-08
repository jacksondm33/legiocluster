process REMOVE_POLY_GS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::python=3.10' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'python-legiocluster:latest' :
        'python-legiocluster:latest' }"

    input:
    tuple val(meta), path(reads)
    val xg

    output:
    tuple val(meta), path("${prefix}_nog_[12].fastq"), emit: nog_reads
    tuple val(meta), path(log_file)                  , emit: log
    path  "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    log_level = "INFO"
    nog_reads = "${prefix}_nog_1.fastq ${prefix}_nog_2.fastq"
    log_file  = "${prefix}.log"

    template 'remove_poly_gs.py'
}
