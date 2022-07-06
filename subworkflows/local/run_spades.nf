include { SPADES              } from '../../modules/local/spades'
include { FILTER_CONTIGS      } from '../../modules/local/filter_contigs'
include { PARSE_SPADES_OUTPUT } from '../../modules/local/parse_spades_output'

workflow RUN_SPADES {
    take:
    reads        // channel: [ val(meta), [ reads ] ]
    max_read_len // channel: [ val(meta), val(max_read_len) ]

    main:
    ch_reports = Channel.empty()
    ch_versions = Channel.empty()

    SPADES (
        reads.join(max_read_len)
    )

    FILTER_CONTIGS (
        SPADES.out.contigs
    )

    PARSE_SPADES_OUTPUT (
        SPADES.out.contigs
    )

    // Collect reports
    ch_reports = ch_reports.concat(PARSE_SPADES_OUTPUT.out.report)

    // Collect versions
    ch_versions = ch_versions.mix(SPADES.out.versions)
    ch_versions = ch_versions.mix(FILTER_CONTIGS.out.versions)
    ch_versions = ch_versions.mix(PARSE_SPADES_OUTPUT.out.versions)

    emit:
    reports = ch_reports
    versions = ch_versions // channel: [ versions.yml ]
}
