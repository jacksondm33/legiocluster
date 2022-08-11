process BWA_MEM {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::bwa=0.7.17" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bwa:0.7.17' :
        'staphb/bwa:0.7.17' }"

    input:
    tuple val(meta), path(reads), path(bwa)

    output:
    tuple val(meta), path(output), emit: sam
    path  "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.ref}"
    output = "${prefix}.sam"
    """
    INDEX=`find -L $bwa -name "*.amb" | sed 's/.amb//'`
    bwa mem \\
        -t $task.cpus \\
        -o $output \\
        $args \\
        \$INDEX \\
        $reads

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """
}
