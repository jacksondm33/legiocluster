include { UNZIP } from '../../modules/nf-core/modules/unzip/main'
include { FASTQC } from '../../modules/local/fastqc'
include { EXTRACT_FASTQC_RESULTS } from '../../modules/local/extract_fastqc_results'
include { CALCULATE_COVERAGE } from '../../modules/local/calculate_coverage'

workflow RUN_FASTQC {
    take:
    raw_reads
    proc_reads
    report

    main:
    ch_versions = Channel.empty()

    FASTQC (
        proc_reads
    )

    UNZIP (
        FASTQC.out.zip.transpose()
    )

    EXTRACT_FASTQC_RESULTS (
        raw_reads.join(UNZIP.out.unzipped_archive.groupTuple(size: 2)).join(report)
    )

    CALCULATE_COVERAGE (
        proc_reads.join(UNZIP.out.unzipped_archive.groupTuple(size: 2)).join(EXTRACT_FASTQC_RESULTS.out.report)
    )

    // Collect versions
    ch_versions = ch_versions.mix(FASTQC.out.versions)
    ch_versions = ch_versions.mix(UNZIP.out.versions)
    ch_versions = ch_versions.mix(EXTRACT_FASTQC_RESULTS.out.versions)
    ch_versions = ch_versions.mix(CALCULATE_COVERAGE.out.versions)

    emit:
    report = CALCULATE_COVERAGE.out.report
    versions = ch_versions // channel: [ versions.yml ]
}
