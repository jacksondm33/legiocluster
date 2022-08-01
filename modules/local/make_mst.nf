process MAKE_MST {
    tag "$meta.ref"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::multiqc=1.12' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.12--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(concat_pairwise_diffs)
    val abr

    output:
    tuple val(meta), path(mst)     , emit: png
    tuple val(meta), path(log_file), emit: log
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.ref}"

    log_level = "INFO"
    mst       = "${prefix}_MST_${abr}.png"
    log_file  = "${prefix}.log"

    template 'make_mst.py'
}
