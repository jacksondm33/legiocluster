include { MASH_SKETCH as MASH_SKETCH_STRAINS } from '../../modules/local/mash_sketch'
include { CONCATENATE                        } from '../../modules/local/concatenate'
include { MASH_SKETCH as MASH_SKETCH_FA      } from '../../modules/local/mash_sketch'
include { MASH_DIST                          } from '../../modules/local/mash_dist'
include { PARSE_MASH_OUTPUT                  } from '../../modules/local/parse_mash_output'

workflow RUN_MASH_FA {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    fasta // channel: [ val(meta), fasta ]

    main:
    ch_reports = Channel.empty()
    ch_versions = Channel.empty()

    MASH_SKETCH_STRAINS (
        Channel.fromPath(params.strain_refs).collect().map { [ [:], it ] }
    )

    MASH_SKETCH_FA (
        fasta
    )

    MASH_DIST (
        MASH_SKETCH_FA.out.mash,
        MASH_SKETCH_STRAINS.out.mash.map { it[1] }
    )

    PARSE_MASH_OUTPUT (
        MASH_DIST.out.dist,
        Channel.fromPath('NO_SPECIES').first()
    )

    // Collect reports
    ch_reports = ch_reports.concat(PARSE_MASH_OUTPUT.out.report)

    // Collect versions
    ch_versions = ch_versions.mix(MASH_SKETCH_STRAINS.out.versions)
    ch_versions = ch_versions.mix(MASH_SKETCH_FA.out.versions)
    ch_versions = ch_versions.mix(MASH_DIST.out.versions)
    ch_versions = ch_versions.mix(PARSE_MASH_OUTPUT.out.versions)

    emit:
    references = PARSE_MASH_OUTPUT.out.references.splitCsv()
    reports = ch_reports
    versions = ch_versions // channel: [ versions.yml ]
}
