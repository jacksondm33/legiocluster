process MASH_SKETCH {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::mash=2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mash:2.1' :
        'staphb/mash:2.1' }"

    input:
    tuple val(meta), path(reads)
    val suffix

    output:
    tuple val(meta), path("*.msh"), emit: mash
    tuple val(meta), path("*.log"), emit: log
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output = "${prefix}_${suffix}.msh"
    def add_params =
        suffix == 'ref_RvSp'   ? '-k 16 -s 400'      :
        suffix == 'comb_reads' ? '-m 2 -k 16 -s 400' :
        ''
    """
    mash \\
        sketch \\
        -o $output \\
        $add_params \\
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
