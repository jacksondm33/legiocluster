process MASH_DIST {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::mash=2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mash:2.1' :
        'staphb/mash:2.1' }"

    input:
    tuple val(meta), path(query), path(reference)

    output:
    tuple val(meta), path(output), emit: dist
    path  "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    suffix = task.ext.suffix ?: 'mash_dist'
    output = "${prefix}.${suffix}.tab"
    """
    mash dist \\
        -p $task.cpus \\
        $args \\
        $reference \\
        $query \\
        > $output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$(mash --version 2>&1)
    END_VERSIONS
    """
}
