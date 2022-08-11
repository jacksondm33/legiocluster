process SAMTOOLS_SORT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::samtools=1.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.9' :
        'staphb/samtools:1.9' }"

    input:
    tuple val(meta), path(sam)

    output:
    tuple val(meta), path(output), emit: bam
    path  "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.ref}"
    output = "${prefix}.bam"
    """
    samtools sort \\
        -@ $task.cpus \\
        -o $output \\
        -T $prefix \\
        $args \\
        $sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
