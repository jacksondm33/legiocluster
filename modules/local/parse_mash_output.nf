process PARSE_MASH_OUTPUT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::python=3.10' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'python-legiocluster:latest' :
        'python-legiocluster:latest' }"

    input:
    tuple val(meta), path(dist)
    path species
    val genome
    val suffix

    output:
    tuple val(meta), path(fastas)  , emit: fastas
    tuple val(meta), path(report)  , emit: report
    tuple val(meta), path(log_file), emit: log
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    log_level = "INFO"
    fastas    = "${prefix}_fastas.csv"
    report    = "${prefix}_mash_${suffix}_report.txt"
    log_file  = "${prefix}.log"

    template 'parse_mash_output.py'
}
