process COUNT_NNN_GAPS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::matplotlib=3.1.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/matplotlib:3.1.2' :
        'quay.io/biocontainers/matplotlib:3.1.2' }"

    input:
    tuple val(meta), path(depth), val(percent_mapped)
    val min_depth
    val gap_length
    val interval
    val max_no_ns
    val max_no_gaps

    output:
    tuple val(meta), path("*_histo_depths.png"), emit: histo_depths
    tuple val(meta), path("*_plot_depths.png") , emit: plot_depths
    tuple val(meta), path("*_report.txt")      , emit: report
    tuple val(meta), path("*.csv")             , emit: csv
    tuple val(meta), path("*.log")             , emit: log
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/legiocluster/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mapped_threshold_param =
        meta.set_ref != 'NO_FILE' ? '--mapped-threshold 0' :
        meta.make_ref == 'true' ? '--mapped-threshold 100' :
        "--mapped-threshold $params.mapped_threshold"
    """
    count_nnn_gaps.py \\
        --depth-file $depth \\
        --histo-depths-file ${prefix}_histo_depths.png \\
        --plot-depths-file ${prefix}_plot_depths.png \\
        --output-file ${prefix}.csv \\
        --report-file ${prefix}_report.txt \\
        --percent-mapped $percent_mapped \\
        --min-depth $min_depth \\
        --gap-length $gap_length \\
        --interval $interval \\
        --max-no-ns $max_no_ns \\
        --max-no-gaps $max_no_gaps \\
        $mapped_threshold_param \\
        $args \\
        > ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
