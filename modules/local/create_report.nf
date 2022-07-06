process CREATE_REPORT {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    tuple val(meta), path("report_*.txt"), path(reads)

    output:
    tuple val(meta), path("report.txt"), emit: report

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def output = "report.txt"
    """
    cat <<EOF > $output
    REPORT

    LegioCluster version:\t\t${workflow.manifest.version}
    Date submitted:\t\t${workflow.start.format(java.time.format.DateTimeFormatter.ISO_LOCAL_DATE)}
    Submitted by:\t\t${workflow.userName}
    Isolate name:\t\t${meta.id}
    Species:\t\t${params.sp_abbr}
    Forward reads:\t\t${reads[0]}
    Reverse reads:\t\t${reads[1]}
    Metadata:\t\t
    Folder name:\t\t${workflow.launchDir}

    EOF
    cat \\
        $args \\
        report_*.txt \\
        >> $output
    """
}
