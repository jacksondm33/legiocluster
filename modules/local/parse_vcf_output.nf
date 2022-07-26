process PARSE_VCF_OUTPUT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::matplotlib=3.1.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/matplotlib:3.1.2' :
        'quay.io/biocontainers/matplotlib:3.1.2' }"

    input:
    tuple val(meta), path(vcf), path(fasta), val(snp_threshold)

    output:
    tuple val(meta), path("*_report.txt"), emit: report
    tuple val(meta), path("*.csv")       , emit: csv
    tuple val(meta), path("*.log")       , emit: log
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/legiocluster/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    parse_vcf_output.py \\
        --isolate $meta.id \\
        --vcf-file $vcf \\
        --reference-file $fasta \\
        --mutation-dist-file ${prefix}_mutation_dist.png \\
        --output-file ${prefix}.csv \\
        --report-file ${prefix}_report.txt \\
        --snp-threshold $snp_threshold \\
        $args \\
        > ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
