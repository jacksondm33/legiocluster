process CONVERT_REPORTS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::python=3.10' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'python-legiocluster:latest' :
        'python-legiocluster:latest' }"

    input:
    tuple val(meta), path(reports)

    output:
    tuple val(meta), path("${prefix}.*.html"), emit: html
    tuple val(meta), path(log_file)          , emit: log
    path  "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    log_level = "INFO"
    log_file  = "${prefix}.log"

    template 'convert_reports.py'
}
