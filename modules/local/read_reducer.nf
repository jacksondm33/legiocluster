process READ_REDUCER {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(reads)
    val both_surviving

    output:
    tuple val(meta), path("*_reduced_*.fastq"), emit: reduced_reads
    tuple val(meta), path("*.log")            , optional:true, emit: log
    path "versions.yml"                       , optional:true, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/legiocluster/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output = "${prefix}_reduced_1.fastq ${prefix}_reduced_2.fastq"
    if (both_surviving < params.min_reads) {
        error "Not enough reads surviving after Trimmomatic."
    }
    if (both_surviving > params.read_cutoff) {
        """
        read_reducer.py \\
            $reads \\
            $output \\
            $args \\
            > ${prefix}.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
        END_VERSIONS
        """
    } else {
        """
        cp ${reads[0]} ${output.split()[0]}
        cp ${reads[1]} ${output.split()[1]}
        """
    }
}
