process REMOVE_POLY_GS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::multiqc=1.12' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.12--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)
    val xg

    output:
    tuple val(meta), path("${prefix}_nog_[12].fastq"), emit: nog_reads
    tuple val(meta), path(log_file)                  , emit: log
    path  "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    log_level = "INFO"
    nog_reads = "${prefix}_nog_1.fastq ${prefix}_nog_2.fastq"
    log_file  = "${prefix}.log"

    template 'remove_poly_gs.py'
}
