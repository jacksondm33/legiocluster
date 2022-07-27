process COMPARE_SNPS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::multiqc=1.12' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.12--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(vcfs), path(csv)

    output:
    tuple val(meta), path("*_SNP_matrix.csv")          , emit: snp_matrix
    tuple val(meta), path("*_ME_matrix.csv")           , emit: me_matrix
    tuple val(meta), path("*_concat_pairwise_snps.csv"), emit: concat_pairwise_snps
    tuple val(meta), path("*_concat_pairwise_mes.csv") , emit: concat_pairwise_mes
    tuple val(meta), path("*.log")                     , emit: log
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/legiocluster/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    log_level            = "INFO"
    log_file             = "${prefix}.log"
    snp_matrix           = "${prefix}_SNP_matrix.csv"
    me_matrix            = "${prefix}_ME_matrix.csv"
    concat_pairwise_snps = "${prefix}_concat_pairwise_snps.csv"
    concat_pairwise_mes  = "${prefix}_concat_pairwise_mes.csv"

    template 'compare_snps.py'
}
