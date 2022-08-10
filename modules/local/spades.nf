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
    tuple val(meta), path('*_contigs.fa'), emit: contigs
    tuple val(meta), path('*.log')       , emit: log
    path  "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def maxmem = task.memory.toGiga()
    def k_param = max_read_len > 175 ? '-k 21,33,55,77,99,127' : '-k 21,33,55,77'
    """
    spades.py \\
        --threads $task.cpus \\
        --memory $maxmem \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        -o ./ \\
        $k_param \\
        $args

    mv spades.log ${prefix}.log
    mv contigs.fasta ${prefix}_SPAdes_contigs.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spades: \$(spades.py --version 2>&1 | sed 's/^.*SPAdes v//; s/ .*\$//')
    END_VERSIONS
    """
}
