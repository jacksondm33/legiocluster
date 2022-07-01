process EXTRACT_FASTQC_RESULTS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(raw_reads), path(fastqc_results), path(report)

    output:
    tuple val(meta), path("*.log")       , emit: log
    tuple val(meta), path(report)        , emit: report
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/legiocluster/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    extract_fastqc_results.py \\
        --fastqc-results ${fastqc_results[0]} \\
        --reads-file ${raw_reads[0]} \\
        --report $report \\
        $args \\
        > ${prefix}.log

    extract_fastqc_results.py \\
        --fastqc-results ${fastqc_results[1]} \\
        --reads-file ${raw_reads[1]} \\
        --report $report \\
        $args \\
        >> ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
