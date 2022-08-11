process PARSE_BWA_OUTPUT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::python=3.10' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'python-legiocluster:latest' :
        'python-legiocluster:latest' }"

    input:
    tuple val(meta), path(fasta), path(sam), path(flagstat), path(idxstats), val(mapped_threshold)

    output:
    tuple val(meta), path(output)  , emit: csv
    tuple val(meta), path(report)  , emit: report
    tuple val(meta), path(log_file), emit: log
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.${meta.ref}"

    log_level = "INFO"
    output    = "${prefix}.csv"
    report    = "${prefix}.bwa_report.txt"
    log_file  = "${prefix}.log"

    template 'parse_bwa_output.py'
}
