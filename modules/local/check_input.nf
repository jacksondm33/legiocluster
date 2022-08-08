process CHECK_INPUT {
    tag "$input"

    conda (params.enable_conda ? 'bioconda::python=3.10' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'python-legiocluster:latest' :
        'python-legiocluster:latest' }"

    input:
    path input
    val references

    output:
    path input_valid   , emit: csv
    path "versions.yml", emit: versions

    script:
    prefix = task.ext.prefix ?: references ? "references" : "samples"

    log_level         = "INFO"
    input_valid       = "${prefix}_valid.csv"
    log_file          = "${prefix}.log"

    template 'check_input.py'
}
