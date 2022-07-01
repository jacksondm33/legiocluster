process CONCATENATE_READS {
    tag "$archive"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("comb_reads.fastq"), emit: comb_reads

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def output = "comb_reads.fastq"
    """
    cat \\
        $reads \\
        $args \\
        > $output
    """
}
