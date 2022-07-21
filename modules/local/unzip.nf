process UNZIP {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::p7zip=15.09" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/p7zip:15.09--h2d50403_4' :
        'quay.io/biocontainers/p7zip:15.09--h2d50403_4' }"

    input:
    tuple val(meta), path(archives)

    output:
    tuple val(meta), path("*_unzipped", type: 'dir'), emit: unzipped_archives
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output = "${prefix}_1_unzipped ${prefix}_2_unzipped"
    """
    7za \\
        e \\
        -o${output.split()[0]}/ \\
        $args \\
        ${archives[0]}

    7za \\
        e \\
        -o${output.split()[1]}/ \\
        $args \\
        ${archives[1]}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        7za: \$(echo \$(7za --help) | sed 's/.*p7zip Version //; s/(.*//')
    END_VERSIONS
    """
}
