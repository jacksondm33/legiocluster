include { FASTQC as FASTQC_MODULE } from '../../modules/local/fastqc'
include { EXTRACT_FASTQC_RESULTS  } from '../../modules/local/extract_fastqc_results'
include { CALCULATE_COVERAGE      } from '../../modules/local/calculate_coverage'

workflow FASTQC {
    take:
    raw_reads  // channel: [ val(meta), [ raw_reads ] ]
    proc_reads // channel: [ val(meta), [ proc_reads ] ]

    main:
    ch_reports = Channel.empty()
    ch_versions = Channel.empty()

    FASTQC_MODULE (
        proc_reads
    )

    EXTRACT_FASTQC_RESULTS (
        FASTQC_MODULE.out.results.join(raw_reads)
    )

    CALCULATE_COVERAGE (
        proc_reads
            .join(FASTQC_MODULE.out.results)
            .map {
                meta, reads, fastqc_results ->
                [ meta, reads[1], fastqc_results[1] ]
            },
        params.med_genome_len
    )

    // Collect reports
    ch_reports = ch_reports.concat(EXTRACT_FASTQC_RESULTS.out.report)
    ch_reports = ch_reports.concat(CALCULATE_COVERAGE.out.report)

    // Collect versions
    ch_versions = ch_versions.mix(FASTQC_MODULE.out.versions)
    ch_versions = ch_versions.mix(EXTRACT_FASTQC_RESULTS.out.versions)
    ch_versions = ch_versions.mix(CALCULATE_COVERAGE.out.versions)

    emit:
    results = FASTQC_MODULE.out.results
    reports = ch_reports
    versions = ch_versions // channel: [ versions.yml ]
}
