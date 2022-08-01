process MASH_SKETCH {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::mash=2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mash:2.1' :
        'staphb/mash:2.1' }"

    input:
    tuple val(meta), path(reads)
    val use_m
    val use_k_s

    output:
    tuple val(meta), path("*.msh"), emit: mash
    tuple val(meta), path("*.log"), emit: log
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output = "${prefix}.msh"
    def m_param = use_m ? '-m 2' : ''
    def k_s_param = use_k_s ? '-k 16 -s 400' : ''
    """
    mash \\
        sketch \\
        -o $output \\
        $m_param \\
        $k_s_param \\
        $args \\
        $reads \\
        > ${prefix}.log

    mash \\
        info \\
        $output \\
        >> ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$(mash --version 2>&1)
    END_VERSIONS
    """
}
