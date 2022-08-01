process COMPARE_SNPS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::multiqc=1.12' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.12--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0' }"

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
    pairwise_diffs = "${prefix}_pairwise_diffs.csv"
    log_file       = "${prefix}.log"

    template 'compare_snps.py'
}