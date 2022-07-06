include { FASTQC                 } from '../../modules/local/fastqc'
include { UNZIP                  } from '../../modules/local/unzip'
include { EXTRACT_FASTQC_RESULTS } from '../../modules/local/extract_fastqc_results'
include { CALCULATE_COVERAGE     } from '../../modules/local/calculate_coverage'

workflow RUN_FASTQC {
    take:
    raw_reads  // channel: [ val(meta), [ raw_reads ] ]
    proc_reads // channel: [ val(meta), [ proc_reads ] ]

    main:
    ch_reports = Channel.empty()
    ch_versions = Channel.empty()

    FASTQC (
        proc_reads
    )

    UNZIP (
        FASTQC.out.zip
    )

    EXTRACT_FASTQC_RESULTS (
        raw_reads.join(UNZIP.out.unzipped_archives)
    )

    CALCULATE_COVERAGE (
        proc_reads.join(UNZIP.out.unzipped_archives)
    )

    // Collect reports
    ch_reports = ch_reports.concat(EXTRACT_FASTQC_RESULTS.out.report)
    ch_reports = ch_reports.concat(CALCULATE_COVERAGE.out.report)

    // Collect versions
    ch_versions = ch_versions.mix(FASTQC.out.versions)
    ch_versions = ch_versions.mix(UNZIP.out.versions)
    ch_versions = ch_versions.mix(EXTRACT_FASTQC_RESULTS.out.versions)
    ch_versions = ch_versions.mix(CALCULATE_COVERAGE.out.versions)

    emit:
    reports = ch_reports
    versions = ch_versions // channel: [ versions.yml ]
}
