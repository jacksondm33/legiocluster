process FILTER_CONTIGS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::python=3.10' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'python-legiocluster:latest' :
        'python-legiocluster:latest' }"

    input:
    tuple val(meta), path(contigs)
    val min_contig_len
    val min_contig_cov

    output:
    tuple val(meta), path(filtered_contigs), emit: filtered_contigs
    tuple val(meta), path(log_file)        , emit: log
    path  "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    log_level        = "INFO"
    filtered_contigs = "${prefix}.fa"
    log_file         = "${prefix}.log"

    template 'filter_contigs.py'
}
