process MAKE_MST {
    tag "$meta.ref"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::python=3.10' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'python-legiocluster:latest' :
        'python-legiocluster:latest' }"

    input:
    tuple val(meta), path(concat_pairwise_diffs)
    val abr
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

    log_level = "INFO"
    mst       = "${prefix}_MST_${abr}.png"
    report    = "${prefix}_report.txt"
    log_file  = "${prefix}.log"

    template 'make_mst.py'
}
