process MASH_SKETCH {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::mash=2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mash:2.1' :
        'staphb/mash:2.1' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path(output)  , emit: mash
    tuple val(meta), path(log_file), emit: log
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    args = task.ext.args ?: ''
    args2 = task.ext.args2 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    suffix = task.ext.suffix ?: 'mash_sketch'
    output = "${prefix}.${suffix}.msh"
    log_file = "${prefix}.${suffix}.log"
    """
    mash sketch \\
        -o $output \\
        $args \\
        $reads \\
        > $log_file

    mash info \\
        $args2 \\
        $output \\
        >> $log_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$(mash --version 2>&1)
    END_VERSIONS
    """
}
