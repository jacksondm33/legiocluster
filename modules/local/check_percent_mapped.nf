process CHECK_PERCENT_MAPPED {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::multiqc=1.12' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.12--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0' }"

    input:
    tuple val(meta), val(percent_mapped), val(min_percent_mapped)

    output:
    tuple val(meta), path(log_file), emit: log
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    log_level = "INFO"
    log_file  = "${prefix}.log"

    template 'check_percent_mapped.py'
}
