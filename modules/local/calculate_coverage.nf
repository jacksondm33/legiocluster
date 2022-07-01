process CALCULATE_COVERAGE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(proc_reads), path(fastqc_results), path(report)

    output:
    tuple val(meta), path("*.log")       , emit: log
    tuple val(meta), path(report)        , emit: report
    tuple val(meta), env(COVERAGE)       , emit: coverage
    tuple val(meta), env(PERC_GE_Q30)    , emit: perc_ge_q30
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/legiocluster/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    calculate_coverage.py \\
        --fastqc-results ${fastqc_results[1]} \\
        --med-genome-len $params.med_genome_len \\
        --reads-file ${proc_reads[1]} \\
        --report-file $report \\
        --summary-file ${prefix}_summary.txt \\
        $args \\
        > ${prefix}.log

    COVERAGE=\$(cat ${prefix}_summary.txt | awk '{print \$1}')
    PERC_GE_Q30=\$(cat ${prefix}_summary.txt | awk '{print \$2}')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
