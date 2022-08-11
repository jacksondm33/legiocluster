process NW_DISPLAY {
    tag "$meta.ref"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::newick_utils=1.6" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/newick_utils:1.6--014d613' :
        'nanozoo/newick_utils:1.6--014d613' }"

    input:
    tuple val(meta), path(parsnp)

    output:
    tuple val(meta), path(output), emit: svg
    path  "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.ref}"
    output = "${prefix}.parsnp_tree_labels.svg"
    VERSION = '1.6'
    """
    nw_display \\
        $args \\
        $parsnp \\
        > $output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        newick_utils: $VERSION
    END_VERSIONS
    """
}
