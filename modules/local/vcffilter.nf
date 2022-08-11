process VCFFILTER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::vcflib=1.0.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vcflib:1.0.1--hd2e4403_1':
        'quay.io/biocontainers/vcflib:1.0.1--hd2e4403_1' }"

    input:
    tuple val(meta), path(vcf), val(dp_max)
    val qual_threshold
    val dp_min
    val qa_threshold
    val ao_dp_ratio

    output:
    tuple val(meta), path(output), emit: vcf
    path  "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    output = "${prefix}.freebayes.vcf"
    VERSION = '1.0.1'
    """
    vcffilter \\
        -f "QUAL > $qual_threshold & DP > $dp_min & DP < $dp_max & QA > $qa_threshold & SAF > 0 & SAR > 0 & AO > $ao_dp_ratio * DP" \\
        $vcf \\
        > $output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcflib: $VERSION
    END_VERSIONS
    """
}
