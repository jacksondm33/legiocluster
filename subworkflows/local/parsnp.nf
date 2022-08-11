include { PARSNP as PARSNP_MODULE } from '../../modules/local/parsnp'
include { PARSE_PARSNP_OUTPUT     } from '../../modules/local/parse_parsnp_output'
include { NW_DISPLAY              } from '../../modules/local/nw_display'
include { REMOVE_INNER_LABELS     } from '../../modules/local/remove_inner_labels'

workflow PARSNP {
    take:
    fasta  // channel: [ meta(ref), fasta      ]
    fastas // channel: [ meta(ref), [ fastas ] ]

    main:
    ch_reports = Channel.empty()
    ch_versions = Channel.empty()

    PARSNP_MODULE (
        fasta.join(fastas)
    )

    PARSE_PARSNP_OUTPUT (
        PARSNP_MODULE.out.parsnp
    )

    NW_DISPLAY (
        PARSNP_MODULE.out.parsnp
    )

    REMOVE_INNER_LABELS (
        NW_DISPLAY.out.svg
    )

    // Collect reports
    ch_reports = ch_reports.concat(PARSE_PARSNP_OUTPUT.out.report)

    // Collect versions
    ch_versions = ch_versions.mix(PARSNP_MODULE.out.versions)
    ch_versions = ch_versions.mix(PARSE_PARSNP_OUTPUT.out.versions)
    ch_versions = ch_versions.mix(NW_DISPLAY.out.versions)
    ch_versions = ch_versions.mix(REMOVE_INNER_LABELS.out.versions)

    emit:
    parsnp = PARSNP_MODULE.out.parsnp
    svg = REMOVE_INNER_LABELS.out.no_inner_labels_svg
    reports = ch_reports
    versions = ch_versions // channel: [ versions.yml ]
}
