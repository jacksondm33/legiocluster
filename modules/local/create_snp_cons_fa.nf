process CREATE_SNP_CONS_FA {
    tag "$meta.ref"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::multiqc=1.12' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.12--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path(snp_cons), emit: snp_cons
    tuple val(meta), path(bases)   , emit: bases
    tuple val(meta), path(log_file), emit: log
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.ref}"

    log_level = "INFO"
    snp_cons  = "${prefix}_SNP_cons.txt"
    bases     = "${prefix}_bases.csv"
    log_file  = "${prefix}.log"

    template 'create_snp_cons_fa.py'
}
