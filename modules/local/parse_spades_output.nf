process PARSE_SPADES_OUTPUT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::python=3.10' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'python-legiocluster:latest' :
        'python-legiocluster:latest' }"

    input:
    tuple val(meta), path(contigs)
    val min_contig_len
    val min_contig_cov
    val max_no_contigs

    output:
    tuple val(meta), path(contig_len_dist)      , emit: contig_len_dist
    tuple val(meta), path(contig_cov_dist)      , emit: contig_cov_dist
    tuple val(meta), path(contig_len_x_cov_dist), emit: contig_len_x_cov_dist
    tuple val(meta), path(contig_ind_len)       , emit: contig_ind_len
    tuple val(meta), path(contig_ind_cov)       , emit: contig_ind_cov
    tuple val(meta), path(report)               , emit: report
    tuple val(meta), path(log_file)             , emit: log
    path  "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    log_level             = "INFO"
    contig_len_dist       = "${prefix}_contig_len_dist.png"
    contig_cov_dist       = "${prefix}_contig_cov_dist.png"
    contig_len_x_cov_dist = "${prefix}_Ampel_dist.png"
    contig_ind_len        = "${prefix}_plot_contig_len.png"
    contig_ind_cov        = "${prefix}_plot_contig_cov.png"
    report                = "${prefix}_report.txt"
    log_file              = "${prefix}.log"

    template 'parse_spades_output.py'
}
