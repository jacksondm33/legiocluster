process FASTQC {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::fastqc=0.11.8" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastqc:0.11.8' :
        'staphb/fastqc:0.11.8' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_fastqc", type: 'dir')                  , emit: results
    tuple val(meta), path("${prefix}_per_base_quality_[12].png")    , emit: per_base_quality
    tuple val(meta), path("${prefix}_per_sequence_quality_[12].png"), emit: per_sequence_quality
    tuple val(meta), path("*.log")                                  , emit: log
    path  "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    fastqc \\
        --threads $task.cpus \\
        $args \\
        $reads \\
        > ${prefix}.log

    cp ${reads[0].baseName}_fastqc/Images/per_base_quality.png ${prefix}_per_base_quality_1.png
    cp ${reads[1].baseName}_fastqc/Images/per_base_quality.png ${prefix}_per_base_quality_2.png
    cp ${reads[0].baseName}_fastqc/Images/per_sequence_quality.png ${prefix}_per_sequence_quality_1.png
    cp ${reads[1].baseName}_fastqc/Images/per_sequence_quality.png ${prefix}_per_sequence_quality_2.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version | sed -e "s/FastQC v//g" )
    END_VERSIONS
    """
}
