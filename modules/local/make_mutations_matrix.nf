process MAKE_MUTATIONS_MATRIX {
    tag "$meta.ref"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::python=3.10' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'python-legiocluster:latest' :
        'python-legiocluster:latest' }"

    input:
    tuple val(meta), path(cluster_pairwise_diffs), path(mutations_matrix)

    output:
    tuple val(meta), path(mutations_matrix, includeInputs: true), emit: mutations_matrix
    tuple val(meta), path(snp_matrix)                           , emit: snp_matrix
    tuple val(meta), path(me_matrix)                            , emit: me_matrix
    tuple val(meta), path(concat_pairwise_snps)                 , emit: concat_pairwise_snps
    tuple val(meta), path(concat_pairwise_mes)                  , emit: concat_pairwise_mes
    tuple val(meta), path(log_file)                             , emit: log
    path  "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.ref}"

    log_level            = "INFO"
    snp_matrix           = "${prefix}_SNP_matrix.csv"
    me_matrix            = "${prefix}_ME_matrix.csv"
    concat_pairwise_snps = "${prefix}_concat_pairwise_snps.csv"
    concat_pairwise_mes  = "${prefix}_concat_pairwise_mes.csv"
    log_file             = "${prefix}.log"

    template 'make_mutations_matrix.py'
}
