include { CONCATENATE       } from '../../modules/local/concatenate'
include { MASH_SKETCH       } from '../../modules/local/mash_sketch'
include { MASH_DIST         } from '../../modules/local/mash_dist'
include { PARSE_MASH_OUTPUT } from '../../modules/local/parse_mash_output'

workflow MASH_FQ {
    take:
    reads // channel: [ meta(id), [ reads ] ]
    mash  // channel: [ meta(id), mash      ]

    main:
    ch_reports = Channel.empty()
    ch_versions = Channel.empty()

    CONCATENATE (
        reads
    )

    MASH_SKETCH (
        CONCATENATE.out.cat,
        'comb_reads'
    )

    MASH_DIST (
        MASH_SKETCH.out.mash.join(mash),
        'RvSp'
    )

    PARSE_MASH_OUTPUT (
        MASH_DIST.out.dist,
        Channel.fromList(
            params.genomes
                .collect {
                    genome, info ->
                    genome + ',' + info.binomial
                })
            .collectFile(name: "species.csv", newLine: true, sort: true)
            .first(),
        params.genome,
        'RvSp'
    )

    // Collect reports
    ch_reports = ch_reports.concat(PARSE_MASH_OUTPUT.out.report)

    // Collect versions
    ch_versions = ch_versions.mix(MASH_SKETCH.out.versions)
    ch_versions = ch_versions.mix(MASH_DIST.out.versions)
    ch_versions = ch_versions.mix(PARSE_MASH_OUTPUT.out.versions)

    emit:
    reports = ch_reports
    versions = ch_versions // channel: [ versions.yml ]
}
