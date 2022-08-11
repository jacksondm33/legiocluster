include { FASTQC as FASTQC_MODULE } from '../../modules/local/fastqc'
include { PARSE_FASTQC_OUTPUT     } from '../../modules/local/parse_fastqc_output'
include { CALCULATE_COVERAGE      } from '../../modules/local/calculate_coverage'

workflow FASTQC {
    take:
    raw_reads  // channel: [ meta(id), [ raw_reads ]  ]
    proc_reads // channel: [ meta(id), [ proc_reads ] ]

    main:
    ch_reports = Channel.empty()
    ch_versions = Channel.empty()

    FASTQC_MODULE (
        proc_reads
    )

    PARSE_FASTQC_OUTPUT (
        FASTQC_MODULE.out.fastqc.join(raw_reads)
    )

    CALCULATE_COVERAGE (
        proc_reads
            .join(FASTQC_MODULE.out.fastqc)
            .map {
                meta, reads, fastqc ->
                [ meta, reads[1], fastqc[1] ]
            },
        params.med_genome_len
    )

    // Collect reports
    ch_reports = ch_reports.concat(PARSE_FASTQC_OUTPUT.out.report)
    ch_reports = ch_reports.concat(CALCULATE_COVERAGE.out.report)

    // Collect versions
    ch_versions = ch_versions.mix(FASTQC_MODULE.out.versions)
    ch_versions = ch_versions.mix(PARSE_FASTQC_OUTPUT.out.versions)
    ch_versions = ch_versions.mix(CALCULATE_COVERAGE.out.versions)

    emit:
    fastqc = FASTQC_MODULE.out.fastqc
    reports = ch_reports
    versions = ch_versions // channel: [ versions.yml ]
}
