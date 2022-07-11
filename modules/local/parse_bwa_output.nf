process PARSE_BWA_OUTPUT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(fasta), path(sam), path(flagstat), path(idxstats)

    output:
    tuple val(meta), path("*_output.csv")        , emit: percent_mapped
    tuple val(meta), path("*_report.txt")        , emit: report
    tuple val(meta), path("*.log")               , emit: log
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/legiocluster/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.ref}"
    def mapped_threshold_param =
        meta.set_ref != 'NO_FILE' ? '--mapped-threshold 0' :
        meta.make_ref == 'true' ? '--mapped-threshold 100' :
        "--mapped-threshold $params.mapped_threshold"
    """
    parse_bwa_output.py \\
        --flagstat-file $flagstat \\
        --idxstats-file $idxstats \\
        --output-file ${prefix}_output.csv \\
        --reference-file $fasta \\
        --report-file ${prefix}_report.txt \\
        --sam-file $sam \\
        $mapped_threshold_param \\
        $args \\
        > ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
