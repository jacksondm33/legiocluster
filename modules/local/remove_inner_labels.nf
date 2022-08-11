process REMOVE_INNER_LABELS {
    tag "$meta.ref"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::python=3.10' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'python-legiocluster:latest' :
        'python-legiocluster:latest' }"

    input:
    tuple val(meta), path(svg)

    output:
    tuple val(meta), path(output)  , emit: no_inner_labels_svg
    tuple val(meta), path(log_file), emit: log
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.ref}"

    log_level = "INFO"
    output    = "${prefix}.parsnp_tree.svg"
    log_file  = "${prefix}.log"

    template 'remove_inner_labels.py'
}
