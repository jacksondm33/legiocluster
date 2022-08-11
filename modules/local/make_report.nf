process MAKE_REPORT {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    tuple val(meta), path(reports), path(reads)
    val genome

    output:
    tuple val(meta), path(output), emit: report

    when:
    task.ext.when == null || task.ext.when

    script:
    args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    output = "${prefix}.report.txt"
    """
    cat <<EOF > $output
    REPORT

    LegioCluster version:\t\t${workflow.manifest.version}
    Date submitted:\t\t${workflow.start.format(java.time.format.DateTimeFormatter.ISO_LOCAL_DATE)}
    Submitted by:\t\t${workflow.userName}
    Isolate name:\t\t${meta.id}
    Species:\t\t${genome}
    Forward reads:\t\t${reads[0]}
    Reverse reads:\t\t${reads[1]}
    Metadata:\t\t
    Folder name:\t\t${workflow.launchDir}

    EOF
    cat \\
        $args \\
        $reports \\
        >> $output
    """
}
