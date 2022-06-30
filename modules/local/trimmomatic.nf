process TRIMMOMATIC {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::trimmomatic=0.39" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trimmomatic:0.39--hdfd78af_2':
        'quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2' }"

    input:
    tuple val(meta), path(reads)
    path "NexteraPE-PE.fa"

    output:
    tuple val(meta), path("*.paired.trim_*.fastq")  , emit: trimmed_reads
    tuple val(meta), path("*.unpaired.trim_*.fastq"), optional:true, emit: unpaired_reads
    tuple val(meta), path("*.log")                  , emit: log
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def trimmed = meta.single_end ? "SE" : "PE"
    def output = meta.single_end ?
        "${prefix}.SE.paired.trim.fastq" // HACK to avoid unpaired and paired in the trimmed_reads output
        : "${prefix}.paired.trim_1.fastq ${prefix}.unpaired.trim_1.fastq ${prefix}.paired.trim_2.fastq ${prefix}.unpaired.trim_2.fastq"
    def qual_trim = params.trimmers
    """
    trimmomatic \\
        $trimmed \\
        -threads $task.cpus \\
        -trimlog ${prefix}.log \\
        $reads \\
        $output \\
        $qual_trim \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trimmomatic: \$(trimmomatic -version)
    END_VERSIONS
    """
}
