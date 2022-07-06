process PARSE_TRIMMOMATIC_LOG {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path("trimlog.txt")

    output:
    tuple val(meta), path("*_report.txt"), emit: report
    tuple val(meta), env(BOTH_SURVIVING) , emit: both_surviving
    tuple val(meta), env(MAX_READ_LEN)   , emit: max_read_len
    tuple val(meta), path("*.log")       , emit: log
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/legiocluster/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    parse_trimmomatic_log.py \\
        --report-file ${prefix}_report.txt \\
        --summary-file ${prefix}_summary.txt \\
        --trimmomatic-log-file trimlog.txt \\
        $args \\
        > ${prefix}.log

    BOTH_SURVIVING=\$(cat ${prefix}_summary.txt | awk '{print \$1}')
    MAX_READ_LEN=\$(cat ${prefix}_summary.txt | awk '{print \$2}')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
