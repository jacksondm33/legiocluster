process MAKE_SNP_CONS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::python=3.10' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'python-legiocluster:latest' :
        'python-legiocluster:latest' }"

    input:
    tuple val(meta), path(fasta), path(mpileup), path(freebayes)

    output:
    tuple val(meta), path(output)  , emit: csv
    tuple val(meta), path(snp_cons), emit: snp_cons
    tuple val(meta), path(log_file), emit: log
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    log_level = "INFO"
    output    = "${prefix}.csv"
    snp_cons  = "${prefix}.SNP_cons.txt"
    log_file  = "${prefix}.log"

    template 'make_snp_cons.py'
}
