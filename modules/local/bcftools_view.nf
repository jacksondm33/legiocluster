process BCFTOOLS_VIEW {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::bcftools=1.9' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.9--h68d8f2e_9' :
        'quay.io/biocontainers/bcftools:1.9--h68d8f2e_9' }"

    input:
    tuple val(meta), path(bcf)

    output:
    tuple val(meta), path(output), emit: vcf
    path  "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.ref}"
    output = "${prefix}.vcf"
    """
    bcftools view \\
        -o $output \\
        $args \\
        $bcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
