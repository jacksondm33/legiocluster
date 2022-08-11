process FREEBAYES {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::freebayes=0.9.21" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/freebayes:0.9.21' :
        'wgspipeline/freebayes:v0.0.1' }"

    input:
    tuple val(meta), path(bam), path(fasta), path(fai)

    output:
    tuple val(meta), path(output), emit: vcf
    path  "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    output = "${prefix}.freebayes_all.vcf"
    """
    freebayes \\
        -f $fasta \\
        $args \\
        $bam \\
        > $output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freebayes: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g' )
    END_VERSIONS
    """
}
