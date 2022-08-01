include { MAKE_MUTATIONS_MATRIX     } from '../../modules/local/make_mutations_matrix'
include { MAKE_MST as MAKE_MST_MES  } from '../../modules/local/make_mst'
include { MAKE_MST as MAKE_MST_SNPS } from '../../modules/local/make_mst'

workflow MAKE_MST {
    take:
    cluster_pairwise_diffs // channel: [ meta(ref), [ cluster_pairwise_diffs ] ]

    main:
    ch_reports = Channel.empty()
    ch_versions = Channel.empty()

    MAKE_MUTATIONS_MATRIX (
        cluster_pairwise_diffs
    )

    MAKE_MST_MES (
        MAKE_MUTATIONS_MATRIX.out.concat_pairwise_mes,
        'ME'
    )

    MAKE_MST_SNPS (
        MAKE_MUTATIONS_MATRIX.out.concat_pairwise_snps,
        'SNP'
    )

    emit:
    reports = ch_reports
    versions = ch_versions // channel: [ versions.yml ]
}