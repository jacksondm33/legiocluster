process MAKE_REFERENCES {
    tag "$references"
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::multiqc=1.12' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.12--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0' }"

    input:
    tuple val(references), val(clusters), val(fastas), val(snp_cons), val(bwas), val(fais), val(mutations_matrices)
    path output

    output:
    tuple val(meta), path(output, includeInputs: true), emit: references
    tuple val(meta), path(log_file)                   , emit: log
    path  "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "references"

    log_level            = "INFO"
    log_file             = "${prefix}.log"

    template 'make_references.py'
}
