process EXTRACT_FASTQC_RESULTS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(reads), path(fastqc_results)

    output:
    tuple val(meta), path("*_per_base_quality_[12].png")    , emit: per_base_quality
    tuple val(meta), path("*_per_sequence_quality_[12].png"), emit: per_sequence_quality
    tuple val(meta), path("*_report.txt")                   , emit: report
    tuple val(meta), path("*.log")                          , emit: log
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/legiocluster/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    extract_fastqc_results.py \\
        --fastqc-results ${fastqc_results[0]} \\
        --reads-file ${reads[0]} \\
        --report ${prefix}_report.txt \\
        $args \\
        > ${prefix}.log

    extract_fastqc_results.py \\
        --fastqc-results ${fastqc_results[1]} \\
        --reads-file ${reads[1]} \\
        --report ${prefix}_report.txt \\
        $args \\
        >> ${prefix}.log

    ln -s ${fastqc_results[0]}/per_base_quality.png ${prefix}_per_base_quality_1.png
    ln -s ${fastqc_results[0]}/per_sequence_quality.png ${prefix}_per_sequence_quality_1.png
    ln -s ${fastqc_results[1]}/per_base_quality.png ${prefix}_per_base_quality_2.png
    ln -s ${fastqc_results[1]}/per_sequence_quality.png ${prefix}_per_sequence_quality_2.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
