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
    tuple val(meta), path("*_references.csv"), emit: references
    tuple val(meta), path("*_report.txt")    , emit: report
    tuple val(meta), path("*.log")           , emit: log
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/legiocluster/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def species_file_param = species.name != 'NO_SPECIES' ? "--species-file $species" : ''
    """
    parse_mash_output.py \\
        --dist-file $dist \\
        --references-file ${prefix}_references.csv \\
        --report-file ${prefix}_report.txt \\
        --sp-abbr $params.sp_abbr \\
        $species_file_param \\
        $args \\
        > ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
