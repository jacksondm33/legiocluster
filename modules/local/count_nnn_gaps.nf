process COUNT_NNN_GAPS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::python=3.10' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'python-legiocluster:latest' :
        'python-legiocluster:latest' }"

    input:
    tuple val(meta), path(depth), val(percent_mapped), val(max_no_ns), val(max_no_gaps), val(mapped_threshold)
    val min_depth
    val gap_length
    val interval

    output:
    tuple val(meta), path(histo_depths), emit: histo_depths
    tuple val(meta), path(plot_depths) , emit: plot_depths
    tuple val(meta), path(output)      , emit: csv
    tuple val(meta), path(report)      , emit: report
    tuple val(meta), path(log_file)    , emit: log
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    log_level    = "INFO"
    histo_depths = "${prefix}.histo_depths.png"
    plot_depths  = "${prefix}.plot_depths.png"
    output       = "${prefix}.csv"
    report       = "${prefix}.nnn_gaps_report.txt"
    log_file     = "${prefix}.log"

    template 'count_nnn_gaps.py'
}
