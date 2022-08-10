process PARSE_PARSNP_OUTPUT {
    tag "$meta.ref"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::python=3.10' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'python-legiocluster:latest' :
        'python-legiocluster:latest' }"

    input:
    tuple val(meta), path(parsnp)

    output:
    tuple val(meta), path(report)  , emit: report
    tuple val(meta), path(log_file), emit: log
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.ref}"

    log_level    = "INFO"
    report       = "${prefix}_parsnp_report.txt"
    log_file     = "${prefix}.log"

    template 'parse_parsnp_output.py'
}
