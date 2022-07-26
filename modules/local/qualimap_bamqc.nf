process QUALIMAP_BAMQC {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::qualimap=2.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/qualimap:2.2.1' :
        'pegi3s/qualimap:2.2.1' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}")     , emit: results
    tuple val(meta), path("*_qualimap.txt"), emit: report
    path  "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def memory = task.memory.toGiga() + "G"
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    qualimap \\
        --java-mem-size=$memory \\
        bamqc \\
        -nt $task.cpus \\
        -bam $bam \\
        -outdir $prefix \\
        $args

    ln -s ${prefix}/genome_results.txt ${prefix}_qualimap.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qualimap: \$(echo \$(qualimap 2>&1) | sed 's/^.*QualiMap v.//; s/Built.*\$//')
    END_VERSIONS
    """
}
