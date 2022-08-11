process SAMTOOLS_IDXSTATS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::samtools=1.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.9' :
        'staphb/samtools:1.9' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path(output), emit: idxstats
    path  "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.ref}"
    output = "${prefix}.idxstats"
    """
    samtools idxstats \\
        --threads $task.cpus \\
        $args \\
        $bam \\
        > $output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
