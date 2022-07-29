include { MAKE_MUTATIONS_MATRIX } from '../../modules/local/make_mutations_matrix'

workflow MAKE_MST {
    take:
    cluster_pairwise_diffs // channel: [ meta(ref), [ cluster_pairwise_diffs ] ]

    main:
    ch_reports = Channel.empty()
    ch_versions = Channel.empty()

    MAKE_MUTATIONS_MATRIX (
        cluster_pairwise_diffs
    )

    // MAKE_MST_MODULE (
    // )

    emit:
    reports = ch_reports
    versions = ch_versions // channel: [ versions.yml ]
}
