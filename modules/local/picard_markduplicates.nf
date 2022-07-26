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
    tuple val(meta), path("*_marked.bam") , emit: bam
    tuple val(meta), path("*_metrics.txt"), emit: metrics
    tuple val(meta), path("*.log")        , emit: log
    path  "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.ref}"
    def maxmem = task.memory.toGiga()
    """
    picard \\
        -Xmx${maxmem}g \\
        MarkDuplicates \\
        $args \\
        I=$bam \\
        O=${prefix}_marked.bam \\
        M=${prefix}_metrics.txt \\
        > ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard MarkDuplicates --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
