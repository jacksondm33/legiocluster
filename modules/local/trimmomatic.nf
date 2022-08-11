process TRIMMOMATIC {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::trimmomatic=0.39" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trimmomatic:0.39':
        'staphb/trimmomatic:0.39' }"

    input:
    tuple val(meta), path(reads)
    path adapters

    output:
    tuple val(meta), path("${prefix}.paired_[12].fastq")  , emit: paired_reads
    tuple val(meta), path("${prefix}.unpaired_[12].fastq"), emit: unpaired_reads
    tuple val(meta), path("${prefix}.log")                , emit: log
    path  "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    args = task.ext.args ?: ''
    args2 = task.ext.args2 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    output = "${prefix}.paired_1.fastq ${prefix}.unpaired_1.fastq ${prefix}.paired_2.fastq ${prefix}.unpaired_2.fastq"
    """
    trimmomatic PE \\
        -threads $task.cpus \\
        -trimlog ${prefix}.log \\
        $args \\
        $reads \\
        $output \\
        $args2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trimmomatic: \$(trimmomatic -version)
    END_VERSIONS
    """
}
