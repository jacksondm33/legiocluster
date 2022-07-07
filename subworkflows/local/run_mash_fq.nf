include { MASH_SKETCH as MASH_SKETCH_SPECIES } from '../../modules/local/mash_sketch'
include { CONCATENATE                        } from '../../modules/local/concatenate'
include { MASH_SKETCH as MASH_SKETCH_FQ      } from '../../modules/local/mash_sketch'
include { MASH_DIST                          } from '../../modules/local/mash_dist'
include { PARSE_MASH_OUTPUT                  } from '../../modules/local/parse_mash_output'

workflow RUN_MASH_FQ {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_reports = Channel.empty()
    ch_versions = Channel.empty()

    MASH_SKETCH_SPECIES (
        Channel.fromPath(params.species_refs).collect().map { [ [:], it ] }
    )

    CONCATENATE (
        reads
    )

    MASH_SKETCH_FQ (
        CONCATENATE.out.cat
    )

    MASH_DIST (
        MASH_SKETCH_FQ.out.mash,
        MASH_SKETCH_SPECIES.out.mash.map { it[1] }
    )

    PARSE_MASH_OUTPUT (
        MASH_DIST.out.dist,
        Channel.of(params.species.collect { sp_abbr, info -> sp_abbr + ',' + info[0] }).flatten().collectFile(name: "species.csv", newLine: true, sort: true)
    )

    // Collect reports
    ch_reports = ch_reports.concat(PARSE_MASH_OUTPUT.out.report)

    // Collect versions
    ch_versions = ch_versions.mix(MASH_SKETCH_SPECIES.out.versions)
    ch_versions = ch_versions.mix(MASH_SKETCH_FQ.out.versions)
    ch_versions = ch_versions.mix(MASH_DIST.out.versions)
    ch_versions = ch_versions.mix(PARSE_MASH_OUTPUT.out.versions)

    emit:
    reports = ch_reports
    versions = ch_versions // channel: [ versions.yml ]
}
