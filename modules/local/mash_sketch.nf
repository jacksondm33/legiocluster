process MASH_SKETCH {
    tag "$meta.id"
    label 'process_medium'
    conda (params.enable_conda ? "bioconda::mash=2.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mash:2.3--he348c14_1' :
        'quay.io/biocontainers/mash:2.3--he348c14_1' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.msh"), emit: mash
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def param = '-m 2 -k 16 -s 400'
    """
    mash \\
        sketch \\
        -p $task.cpus \\
        -o ${prefix}.msh \\
        $param \\
        $args \\
        $reads \\
        > ${prefix}.log

    mash \\
        info \\
        ${prefix}.msh \\
        >> ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$(mash --version 2>&1)
    END_VERSIONS
    """
}
