process PARSE_MASH_OUTPUT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(dist)
    path species

    output:
    tuple val(meta), path("*_report.txt"), emit: report
    tuple val(meta), env(MASH_SPECIES)   , emit: mash_species
    tuple val(meta), env(PASSED_QC)      , emit: passed_qc
    tuple val(meta), path("*.log")       , emit: log
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/legiocluster/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    parse_mash_output.py \\
        --dist-file $dist \\
        --report-file ${prefix}_report.txt \\
        --species-file $species \\
        --summary-file ${prefix}_summary.txt \\
        $args \\
        > ${prefix}.log

    MASH_SPECIES=\$(cat ${prefix}_summary.txt | awk '{print \$1}')
    PASSED_QC=\$(cat ${prefix}_summary.txt | awk '{print \$2}')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
