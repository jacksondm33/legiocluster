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
    tuple val(meta), path(output, type: 'dir')     , emit: qualimap
    tuple val(meta), path("${prefix}.qualimap.txt"), emit: report
    path  "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    output = "${prefix}"
    """
    qualimap bamqc \\
        -nt $task.cpus \\
        -bam $bam \\
        -outdir $output \\
        $args

    cp ${output}/genome_results.txt ${prefix}.qualimap.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qualimap: \$(echo \$(qualimap 2>&1) | sed 's/^.*QualiMap v.//; s/Built.*\$//')
    END_VERSIONS
    """
}
