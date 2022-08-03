process CHECK_REF_QUAL {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::multiqc=1.12' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.12--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(contigs)
    val med_genome_len

    output:
    tuple val(meta), path(output)  , emit: csv
    tuple val(meta), path(report)  , emit: report
    tuple val(meta), path(log_file), emit: log
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    log_level = "INFO"
    output    = "${prefix}.csv"
    report    = "${prefix}_report.txt"
    log_file  = "${prefix}.log"

    template 'check_ref_qual.py'
}
