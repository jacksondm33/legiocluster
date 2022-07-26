process CLEANUP_FREEBAYES {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::multiqc=1.12' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.12--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta), path(mpileup), path(vcf)
    val diagnostic_mode

    output:
    tuple val(meta), path("*_SNP_cons.txt"), emit: snp_cons
    tuple val(meta), path("*.log")         , emit: log
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/legiocluster/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    log_level = "INFO"
    log_file = "${prefix}.log"
    snp_cons = "${prefix}_SNP_cons.txt"
    output = "${prefix}.csv"

    template 'cleanup_freebayes.py'
}
