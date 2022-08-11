process PARSE_KRAKEN_OUTPUT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::python=3.10' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'python-legiocluster:latest' :
        'python-legiocluster:latest' }"

    input:
    tuple val(meta), path(kraken), path(contigs)
    val genus

    output:
    tuple val(meta), path(good_contigs), emit: good_contigs
    tuple val(meta), path(bad_contigs) , emit: bad_contigs
    tuple val(meta), path(report)      , emit: report
    tuple val(meta), path(log_file)    , emit: log
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    log_level    = "INFO"
    good_contigs = "${prefix}_cc.fasta"
    bad_contigs  = "${prefix}.wrong_genus_contigs.fasta"
    report       = "${prefix}.kraken_report.txt"
    log_file     = "${prefix}.log"

    template 'parse_kraken_output.py'
}
