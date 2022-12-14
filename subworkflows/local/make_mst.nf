include { MAKE_MUTATIONS_MATRIX     } from '../../modules/local/make_mutations_matrix'
include { MAKE_MST as MAKE_MST_ME   } from '../../modules/local/make_mst'
include { MAKE_MST as MAKE_MST_SNP  } from '../../modules/local/make_mst'

workflow MAKE_MST {
    take:
    cluster_pairwise_diffs // channel: [ meta(ref), [ cluster_pairwise_diffs ] ]
    mutations_matrix       // channel: [ meta(ref), mutations_matrix           ]

    main:
    ch_reports = Channel.empty()
    ch_versions = Channel.empty()

    MAKE_MUTATIONS_MATRIX (
        cluster_pairwise_diffs.join(mutations_matrix)
    )

    MAKE_MST_ME (
        MAKE_MUTATIONS_MATRIX.out.concat_pairwise_mes,
        params.genome
    )

    MAKE_MST_SNP (
        MAKE_MUTATIONS_MATRIX.out.concat_pairwise_snps,
        params.genome
    )

    // Collect reports
    ch_reports = ch_reports.concat(MAKE_MST_ME.out.report)
    ch_reports = ch_reports.concat(MAKE_MST_SNP.out.report)

    // Collect versions
    ch_versions = ch_versions.mix(MAKE_MUTATIONS_MATRIX.out.versions)
    ch_versions = ch_versions.mix(MAKE_MST_ME.out.versions)
    ch_versions = ch_versions.mix(MAKE_MST_SNP.out.versions)

    emit:
    mutations_matrix = MAKE_MUTATIONS_MATRIX.out.mutations_matrix
    reports = ch_reports
    versions = ch_versions // channel: [ versions.yml ]
}
