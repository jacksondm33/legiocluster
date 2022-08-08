include { CONCATENATE                           } from '../../modules/local/concatenate'
include { MASH_SKETCH                           } from '../../modules/local/mash_sketch'
include { MASH_SKETCH as MASH_SKETCH_REFERENCES } from '../../modules/local/mash_sketch'
include { MASH_DIST                             } from '../../modules/local/mash_dist'
include { PARSE_MASH_OUTPUT                     } from '../../modules/local/parse_mash_output'

workflow MASH_FA {
    take:
    fasta  // channel: [ meta(id), fasta      ]
    fastas // channel: [ meta(id), [ fastas ] ]

    main:
    ch_reports = Channel.empty()
    ch_versions = Channel.empty()

    MASH_SKETCH (
        fasta,
        'query_FAvNCBI'
    )

    MASH_SKETCH_REFERENCES (
        fastas,
        'ref_FAvNCBI'
    )

    MASH_DIST (
        MASH_SKETCH.out.mash.join(MASH_SKETCH_REFERENCES.out.mash),
        'FAvNCBI'
    )

    PARSE_MASH_OUTPUT (
        MASH_DIST.out.dist,
        Channel.fromPath('NO_FILE').first(),
        params.genome
    )

    // Collect reports
    ch_reports = ch_reports.concat(PARSE_MASH_OUTPUT.out.report)

    // Collect versions
    ch_versions = ch_versions.mix(MASH_SKETCH.out.versions)
    ch_versions = ch_versions.mix(MASH_SKETCH_REFERENCES.out.versions)
    ch_versions = ch_versions.mix(MASH_DIST.out.versions)
    ch_versions = ch_versions.mix(PARSE_MASH_OUTPUT.out.versions)

    emit:
    fastas = PARSE_MASH_OUTPUT.out.fastas.map { meta, csv -> [ meta ] + csv.splitCsv().collect() }
    reports = ch_reports
    versions = ch_versions // channel: [ versions.yml ]
}
