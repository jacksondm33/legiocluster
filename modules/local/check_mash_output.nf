process CHECK_MASH_OUTPUT {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    tuple val(meta), val(mash_species), val(passed_qc)

    output:

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/legiocluster/bin/
    if (passed_qc != "True") {
        error "The reads did not pass the Mash QC check."
    }
    if (mash_species != params.sp_abbr) {
        error "The reads did not pass the Mash species check."
    }
    """
    """
}
