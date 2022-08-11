process PICARD_MARKDUPLICATES {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::picard=2.18.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.18.2--py36_0' :
        'quay.io/biocontainers/picard:2.18.2--py36_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.marked.bam") , emit: marked_bam
    tuple val(meta), path("${prefix}.metrics.txt"), emit: metrics
    path  "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.ref}"
    """
    picard MarkDuplicates \\
        $args \\
        I=$bam \\
        O=${prefix}.marked.bam \\
        M=${prefix}.metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard MarkDuplicates --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
