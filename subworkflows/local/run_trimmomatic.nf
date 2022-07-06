include { GUNZIP                } from '../../modules/local/gunzip'
include { REMOVE_POLY_GS        } from '../../modules/local/remove_poly_gs'
include { TRIMMOMATIC           } from '../../modules/local/trimmomatic'
include { PARSE_TRIMMOMATIC_LOG } from '../../modules/local/parse_trimmomatic_log'
include { REDUCE_READS          } from '../../modules/local/reduce_reads'

workflow RUN_TRIMMOMATIC {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_reports = Channel.empty()
    ch_versions = Channel.empty()

    GUNZIP (
        reads
    )

    REMOVE_POLY_GS (
        GUNZIP.out.unzipped_reads
    )

    TRIMMOMATIC (
        REMOVE_POLY_GS.out.nog_reads,
        Channel.fromPath(params.illuminaclip).first()
    )

    PARSE_TRIMMOMATIC_LOG (
        TRIMMOMATIC.out.log
    )

    REDUCE_READS (
        TRIMMOMATIC.out.trimmed_reads.join(PARSE_TRIMMOMATIC_LOG.out.both_surviving)
    )

    // Collect reports
    ch_reports = ch_reports.concat(PARSE_TRIMMOMATIC_LOG.out.report)

    // Collect versions
    ch_versions = ch_versions.mix(GUNZIP.out.versions)
    ch_versions = ch_versions.mix(REMOVE_POLY_GS.out.versions)
    ch_versions = ch_versions.mix(TRIMMOMATIC.out.versions)
    ch_versions = ch_versions.mix(PARSE_TRIMMOMATIC_LOG.out.versions)
    ch_versions = ch_versions.mix(REDUCE_READS.out.versions)

    emit:
    reads = REDUCE_READS.out.reduced_reads
    max_read_len = PARSE_TRIMMOMATIC_LOG.out.max_read_len
    reports = ch_reports
    versions = ch_versions // channel: [ versions.yml ]
}
