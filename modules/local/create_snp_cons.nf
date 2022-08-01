process CREATE_SNP_CONS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::multiqc=1.12' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.12--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(mpileup), path(vcf), path(bases)
    val diagnostic_mode

    output:
    tuple val(meta), path(csv)     , emit: csv
    tuple val(meta), path(snp_cons), emit: snp_cons
    tuple val(meta), path(log_file), emit: log
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    log_level = "INFO"
    csv       = "${prefix}.csv"
    snp_cons  = "${prefix}_SNP_cons.txt"
    log_file  = "${prefix}.log"

    template 'create_snp_cons.py'
}