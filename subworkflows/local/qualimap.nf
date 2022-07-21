include { QUALIMAP_BAMQC        } from '../../modules/local/qualimap_bamqc'
include { PARSE_QUALIMAP_OUTPUT } from '../../modules/local/parse_qualimap_output'

workflow QUALIMAP {
    take:
    bam  // channel: [ val(meta), bam ]

    main:
    ch_reports = Channel.empty()
    ch_versions = Channel.empty()

    QUALIMAP_BAMQC (
        bam
    )

    PARSE_QUALIMAP_OUTPUT (
        QUALIMAP_BAMQC.out.report
    )

    // Collect reports
    ch_reports = ch_reports.concat(PARSE_QUALIMAP_OUTPUT.out.report)

    // Collect versions
    ch_versions = ch_versions.mix(QUALIMAP_BAMQC.out.versions)
    ch_versions = ch_versions.mix(PARSE_QUALIMAP_OUTPUT.out.versions)

    emit:
    reports = ch_reports
    versions = ch_versions // channel: [ versions.yml ]
}
