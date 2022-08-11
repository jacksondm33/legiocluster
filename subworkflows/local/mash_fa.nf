include { MASH_SKETCH as MASH_SKETCH_QUERY_FA       } from '../../modules/local/mash_sketch'
include { MASH_SKETCH as MASH_SKETCH_REF_FA         } from '../../modules/local/mash_sketch'
include { MASH_DIST as MASH_DIST_FA                 } from '../../modules/local/mash_dist'
include { PARSE_MASH_OUTPUT as PARSE_MASH_OUTPUT_FA } from '../../modules/local/parse_mash_output'

workflow MASH_FA {
    take:
    fasta  // channel: [ meta(id), fasta      ]
    fastas // channel: [ meta(id), [ fastas ] ]

    main:
    ch_reports = Channel.empty()
    ch_versions = Channel.empty()

    MASH_SKETCH_QUERY_FA (
        fasta
    )

    MASH_SKETCH_REF_FA (
        fastas
    )

    MASH_DIST_FA (
        MASH_SKETCH_QUERY_FA.out.mash.join(MASH_SKETCH_REF_FA.out.mash)
    )

    PARSE_MASH_OUTPUT_FA (
        MASH_DIST_FA.out.dist,
        Channel.fromPath('NO_FILE').first(),
        params.genome
    )

    PARSE_MASH_OUTPUT_FA.out.fastas
        .map {
            meta, csv ->
            [ meta ] + csv.splitCsv().collect()
        }
        .set { ch_fastas }

    // Collect reports
    ch_reports = ch_reports.concat(PARSE_MASH_OUTPUT_FA.out.report)

    // Collect versions
    ch_versions = ch_versions.mix(MASH_SKETCH_QUERY_FA.out.versions)
    ch_versions = ch_versions.mix(MASH_SKETCH_REF_FA.out.versions)
    ch_versions = ch_versions.mix(MASH_DIST_FA.out.versions)
    ch_versions = ch_versions.mix(PARSE_MASH_OUTPUT_FA.out.versions)

    emit:
    fastas = ch_fastas
    reports = ch_reports
    versions = ch_versions // channel: [ versions.yml ]
}
