process CREATE_SNP_CONS_FA {
    tag "$meta.ref"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::python=3.10' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'python-legiocluster:latest' :
        'python-legiocluster:latest' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path(snp_cons), emit: snp_cons
    tuple val(meta), path(log_file), emit: log
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.ref}"

    log_level = "INFO"
    snp_cons  = "${prefix}_SNP_cons.txt"
    log_file  = "${prefix}.log"

    template 'create_snp_cons_fa.py'
}
