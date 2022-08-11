process CALCULATE_COVERAGE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::python=3.10' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'python-legiocluster:latest' :
        'python-legiocluster:latest' }"

    input:
    tuple val(meta), path(reads), path(fastqc)
    val med_genome_len

    output:
    tuple val(meta), path(report)  , emit: report
    tuple val(meta), path(log_file), emit: log
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    log_level = "INFO"
    report    = "${prefix}.coverage_report.txt"
    log_file  = "${prefix}.log"

    template 'calculate_coverage.py'
}
