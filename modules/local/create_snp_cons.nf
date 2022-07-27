process CREATE_SNP_CONS {
    tag "$meta.ref"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::multiqc=1.12' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.12--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta), path(vcfs)

    output:
    tuple val(meta), path("*_bases.csv"), emit: bases
    tuple val(meta), path("*.log")      , emit: log
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/legiocluster/bin/
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.ref}"

    log_level = "INFO"
    log_file  = "${prefix}.log"
    bases     = "${prefix}_bases.csv"

    template 'create_snp_cons.py'
}
