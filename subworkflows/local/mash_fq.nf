include { CONCATENATE                               } from '../../modules/local/concatenate'
include { MASH_SKETCH as MASH_SKETCH_FQ             } from '../../modules/local/mash_sketch'
include { MASH_DIST as MASH_DIST_FQ                 } from '../../modules/local/mash_dist'
include { PARSE_MASH_OUTPUT as PARSE_MASH_OUTPUT_FQ } from '../../modules/local/parse_mash_output'

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

    MASH_SKETCH_FQ (
        CONCATENATE.out.cat
    )

    MASH_DIST_FQ (
        MASH_SKETCH_FQ.out.mash.join(mash)
    )

    PARSE_MASH_OUTPUT_FQ (
        MASH_DIST_FQ.out.dist,
        Channel.fromList(
            params.genomes
                .collect {
                    genome, info ->
                    genome + ',' + info.binomial
                })
            .collectFile(name: "species.csv", newLine: true, sort: true)
            .first(),
        params.genome
    )

    // Collect reports
    ch_reports = ch_reports.concat(PARSE_MASH_OUTPUT_FQ.out.report)

    // Collect versions
    ch_versions = ch_versions.mix(MASH_SKETCH_FQ.out.versions)
    ch_versions = ch_versions.mix(MASH_DIST_FQ.out.versions)
    ch_versions = ch_versions.mix(PARSE_MASH_OUTPUT_FQ.out.versions)

    emit:
    reports = ch_reports
    versions = ch_versions // channel: [ versions.yml ]
}
