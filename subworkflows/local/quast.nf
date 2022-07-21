include { CHECK_PERCENT_MAPPED  } from '../../modules/local/check_percent_mapped'
include { QUAST as QUAST_MODULE } from '../../modules/local/quast'
include { PARSE_QUAST_OUTPUT    } from '../../modules/local/parse_quast_output'
include { COUNT_NNN_GAPS        } from '../../modules/local/count_nnn_gaps'

workflow QUAST {
    take:
    contigs  // channel: [ val(meta), [ raw_reads ] ]
    fasta    // channel: [ val(meta), [ proc_reads ] ]
    depth
    percent_mapped

    main:
    ch_reports = Channel.empty()
    ch_versions = Channel.empty()

    CHECK_PERCENT_MAPPED (
        percent_mapped
    )

    QUAST_MODULE (
        contigs.join(fasta)
    )

    PARSE_QUAST_OUTPUT (
        contigs.join(QUAST_MODULE.out.report)
    )

    COUNT_NNN_GAPS (
        depth.join(percent_mapped),
        params.min_depth,
        params.gap_length,
        params.interval,
        params.max_no_ns,
        params.max_no_gaps
    )

    // Collect reports
    ch_reports = ch_reports.concat(PARSE_QUAST_OUTPUT.out.report)
    ch_reports = ch_reports.concat(COUNT_NNN_GAPS.out.report)

    // Collect versions
    ch_versions = ch_versions.mix(QUAST_MODULE.out.versions)
    ch_versions = ch_versions.mix(PARSE_QUAST_OUTPUT.out.versions)
    ch_versions = ch_versions.mix(COUNT_NNN_GAPS.out.versions)

    emit:
    reports = ch_reports
    versions = ch_versions // channel: [ versions.yml ]
}
