process KRAKEN {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? 'bioconda::kraken=1.1.1' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kraken:1.1.1' :
        'staphb/kraken:1.1.1' }"

    input:
    tuple val(meta), path(fasta)
    path db

    output:
    tuple val(meta), path(output), emit: kraken
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    args = task.ext.args ?: ''
    args2 = task.ext.args2 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    output = "${prefix}.kraken_res.txt"
    """
    kraken \\
        --threads $task.cpus \\
        --db $db \\
        --output ${prefix}.kraken_out.txt \\
        $args \\
        $fasta

    kraken-translate \\
        --db $db \\
        $args2 \\
        ${prefix}.kraken_out.txt \\
        > $output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken: \$(echo \$(kraken --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
    END_VERSIONS
    """
}
