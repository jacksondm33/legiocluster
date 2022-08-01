process PARSE_VCF_OUTPUT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::multiqc=1.12' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.12--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0' }"

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
    mutation_dist = "${prefix}_mutation_dist.png"
    output        = "${prefix}.csv"
    report        = "${prefix}_report.txt"
    log_file      = "${prefix}.log"

    template 'parse_vcf_output.py'
}
