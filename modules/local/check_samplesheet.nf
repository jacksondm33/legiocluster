process CHECK_SAMPLESHEET {
    tag "$samplesheet"

    conda (params.enable_conda ? 'bioconda::multiqc=1.12' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.12--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0' }"

    input:
    path samplesheet

    output:
    path samplesheet_valid, emit: csv
    path "versions.yml"   , emit: versions

    script:
    prefix = task.ext.prefix ?: "samplesheet"

    log_level         = "INFO"
    samplesheet_valid = "${prefix}_valid.csv"
    log_file          = "${prefix}.log"

    template 'check_samplesheet.py'
}
