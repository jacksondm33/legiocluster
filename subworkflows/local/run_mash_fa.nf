include { CONCATENATE                        } from '../../modules/local/concatenate'
include { MASH_SKETCH                        } from '../../modules/local/mash_sketch'
include { MASH_DIST as MASH_DIST_FA          } from '../../modules/local/mash_dist'
include { PARSE_MASH_OUTPUT                  } from '../../modules/local/parse_mash_output'

workflow RUN_MASH_FA {
    take:
    reads // channel: [ meta(id), [ reads ] ]
    fasta // channel: [ meta(id), fasta ]
    mash  // channel: [ mash ]

    main:
    ch_reports = Channel.empty()
    ch_versions = Channel.empty()

    MASH_SKETCH (
        fasta
    )

    MASH_DIST_FA (
        MASH_SKETCH.out.mash,
        mash
    )

    PARSE_MASH_OUTPUT (
        MASH_DIST_FA.out.dist,
        Channel.fromPath('NO_FILE').first()
    )

    // Collect reports
    ch_reports = ch_reports.concat(PARSE_MASH_OUTPUT.out.report)

    // Collect versions
    ch_versions = ch_versions.mix(MASH_SKETCH.out.versions)
    ch_versions = ch_versions.mix(MASH_DIST_FA.out.versions)
    ch_versions = ch_versions.mix(PARSE_MASH_OUTPUT.out.versions)

    emit:
    fastas = PARSE_MASH_OUTPUT.out.fastas.map { meta, fastas -> [ meta, fastas.splitCsv().collect { it[0] } ] }
    reports = ch_reports
    versions = ch_versions // channel: [ versions.yml ]
}
