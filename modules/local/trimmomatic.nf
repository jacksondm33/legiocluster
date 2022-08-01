process TRIMMOMATIC {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::trimmomatic=0.39" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trimmomatic:0.39':
        'staphb/trimmomatic:0.39' }"

    input:
    tuple val(meta), path(reads)
    path "NexteraPE-PE.fa"

    output:
    tuple val(meta), path("*_paired_[12].fastq")  , emit: trimmed_reads
    tuple val(meta), path("*_unpaired_[12].fastq"), emit: unpaired_reads
    tuple val(meta), path("*.log")                , emit: log
    path  "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def trimmed = meta.single_end ? "SE" : "PE"
    def output = "${prefix}_paired_1.fastq ${prefix}_unpaired_1.fastq ${prefix}_paired_2.fastq ${prefix}_unpaired_2.fastq"
    def qual_trim = task.ext.args2 ?: ''
    """
    trimmomatic \\
        $trimmed \\
        -threads $task.cpus \\
        -trimlog ${prefix}.log \\
        $args \\
        $reads \\
        $output \\
        $qual_trim

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trimmomatic: \$(trimmomatic -version)
    END_VERSIONS
    """
}
