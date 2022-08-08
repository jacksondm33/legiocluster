include { QUAST as QUAST_MODULE } from '../../modules/local/quast'
include { PARSE_QUAST_OUTPUT    } from '../../modules/local/parse_quast_output'
include { COUNT_NNN_GAPS        } from '../../modules/local/count_nnn_gaps'

workflow QUAST {
    take:
    contigs  // channel: [ val(meta), [ raw_reads ] ]
    fasta    // channel: [ val(meta), [ proc_reads ] ]
    depth
    percent_mapped
    max_no_ns
    max_no_gaps
    mapped_threshold

    main:
    ch_reports = Channel.empty()
    ch_versions = Channel.empty()

    QUAST_MODULE (
        contigs.join(fasta)
    )

    PARSE_QUAST_OUTPUT (
        contigs.join(QUAST_MODULE.out.report),
        params.contig_threshold
    )

    COUNT_NNN_GAPS (
        depth.join(percent_mapped).join(max_no_ns).join(max_no_gaps).join(mapped_threshold),
        params.min_depth,
        params.gap_length,
        params.interval
    )

    COUNT_NNN_GAPS.out.csv
        .map {
            meta, csv ->
            [ meta ] + csv.splitCsv().collect()
        }
        .multiMap {
            meta, depth_mean, depth_sd ->
            depth_mean: [ meta, depth_mean[0].toFloat() ]
            depth_sd: [ meta, depth_sd[0].toFloat() ]
        }
        .set { ch_output }

    // Collect reports
    ch_reports = ch_reports.concat(PARSE_QUAST_OUTPUT.out.report)
    ch_reports = ch_reports.concat(COUNT_NNN_GAPS.out.report)

    // Collect versions
    ch_versions = ch_versions.mix(QUAST_MODULE.out.versions)
    ch_versions = ch_versions.mix(PARSE_QUAST_OUTPUT.out.versions)
    ch_versions = ch_versions.mix(COUNT_NNN_GAPS.out.versions)

    emit:
    depth_mean = ch_output.depth_mean
    depth_sd = ch_output.depth_sd
    reports = ch_reports
    versions = ch_versions // channel: [ versions.yml ]
}
