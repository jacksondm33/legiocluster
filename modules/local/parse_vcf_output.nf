process PARSE_VCF_OUTPUT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::python=3.10' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'python-legiocluster:latest' :
        'python-legiocluster:latest' }"

    input:
    tuple val(meta), path(vcf), path(fasta), val(snp_threshold)

    output:
    tuple val(meta), path(mutation_dist), emit: mutation_dist
    tuple val(meta), path(output)       , emit: csv
    tuple val(meta), path(report)       , emit: report
    tuple val(meta), path(log_file)     , emit: log
    path  "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    log_level     = "INFO"
    mutation_dist = "${prefix}.mutation_dist.png"
    output        = "${prefix}.csv"
    report        = "${prefix}.vcf_report.txt"
    log_file      = "${prefix}.log"

    template 'parse_vcf_output.py'
}
