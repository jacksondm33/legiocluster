process PARSE_SPADES_OUTPUT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::matplotlib=3.1.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/matplotlib:3.1.2' :
        'quay.io/biocontainers/matplotlib:3.1.2' }"

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("*_report.txt"), emit: report
    tuple val(meta), path("*.png")       , emit: plots
    tuple val(meta), path("*.log")       , emit: log
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/legiocluster/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    parse_spades_output.py \\
        --contigs-file $contigs \\
        --contig-len-dist-file ${prefix}_contig_len_dist.png \\
        --contig-cov-dist-file ${prefix}_contig_cov_dist.png \\
        --contig-len-x-cov-dist-file ${prefix}_Ampel_dist.png \\
        --contig-ind-len-file ${prefix}_plot_contig_len.png \\
        --contig-ind-cov-file ${prefix}_plot_contig_cov.png \\
        --report-file ${prefix}_report.txt \\
        $args \\
        > ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
