process QUAST {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::quast=5.0.2' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/quast:5.0.2' :
        'staphb/quast:5.0.2' }"

    input:
    tuple val(meta), path(contigs), path(fasta)

    output:
    tuple val(meta), path("${prefix}")  , emit: results
    tuple val(meta), path("*_quast.txt"), emit: report
    path  "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    quast.py \\
        -t $task.cpus \\
        -o $prefix \\
        -r $fasta \\
        $args \\
        $contigs

    ln -s ${prefix}/report.txt ${prefix}_quast.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: \$(quast.py --version 2>&1 | sed 's/^.*QUAST v//; s/ .*\$//')
    END_VERSIONS
    """
}
