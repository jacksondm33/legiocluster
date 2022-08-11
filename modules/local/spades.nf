process SPADES {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? 'bioconda::spades=3.12.0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/spades:3.12.0' :
        'staphb/spades:3.12.0' }"

    input:
    tuple val(meta), path(reads), val(max_read_len)

    output:
    tuple val(meta), path("${prefix}.SPAdes_contigs.fa"), emit: contigs
    tuple val(meta), path("${prefix}.log")              , emit: log
    path  "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    maxmem = task.memory.toGiga()
    k_param = max_read_len > 175 ? '-k 21,33,55,77,99,127' : '-k 21,33,55,77'
    """
    spades.py \\
        --threads $task.cpus \\
        --memory $maxmem \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        -o ./ \\
        $k_param \\
        $args

    mv contigs.fasta ${prefix}.SPAdes_contigs.fa
    mv spades.log ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spades: \$(spades.py --version 2>&1 | sed 's/^.*SPAdes v//; s/ .*\$//')
    END_VERSIONS
    """
}
