process MAKE_MST {
    tag "$meta.ref"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::python=3.10' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'python-legiocluster:latest' :
        'python-legiocluster:latest' }"

    input:
    tuple val(meta), path(concat_pairwise_diffs)
    val genome

    output:
    tuple val(meta), path(mst)     , emit: png
    tuple val(meta), path(report)  , emit: report
    tuple val(meta), path(log_file), emit: log
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.ref}"
    suffix = task.ext.suffix ?: 'ME'

    log_level = "INFO"
    mst       = "${prefix}.MST_${suffix}.png"
    report    = "${prefix}.MST_${suffix}_report.txt"
    log_file  = "${prefix}.MST_${suffix}.log"

    template 'make_mst.py'
}
