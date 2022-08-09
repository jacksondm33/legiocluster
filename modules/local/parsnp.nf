process PARSNP {
    tag "$meta.ref"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::parsnp=1.5.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/parsnp:1.5.3' :
        'quay.io/biocontainers/parsnp:1.5.3--he513fc3_0' }"

    input:
    tuple val(meta), path(fasta), path(fastas)

    output:
    tuple val(meta), path(output) , emit: tree
    tuple val(meta), path("*.log"), emit: log
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.ref}"
    output = "${prefix}_parsnp.tree"
    """
    parsnp \\
        -p $task.cpus \\
        -o $prefix \\
        -d $fasta $fastas \\
        -r $fasta \\
        $args \\
        > ${prefix}.log

    mv ${prefix}/parsnp.tree $output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parsnp: \$(echo \$(parsnp --version 2>&1) | sed 's/^.*(parsnp) //; s/ .*\$//')
    END_VERSIONS
    """
}
