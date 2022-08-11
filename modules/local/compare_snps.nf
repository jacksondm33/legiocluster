process COMPARE_SNPS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::python=3.10' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'python-legiocluster:latest' :
        'python-legiocluster:latest' }"

    input:
    tuple val(meta), path(snp_cons), path(cluster_snp_cons)

    output:
    tuple val(meta), path(pairwise_diffs), emit: pairwise_diffs
    tuple val(meta), path(log_file)      , emit: log
    path  "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    log_level      = "INFO"
    pairwise_diffs = "${prefix}.pairwise_diffs.csv"
    log_file       = "${prefix}.log"

    template 'compare_snps.py'
}
