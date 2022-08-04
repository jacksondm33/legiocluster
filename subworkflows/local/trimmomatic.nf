include { GUNZIP                            } from '../../modules/local/gunzip'
include { REMOVE_POLY_GS                    } from '../../modules/local/remove_poly_gs'
include { TRIMMOMATIC as TRIMMOMATIC_MODULE } from '../../modules/local/trimmomatic'
include { PARSE_TRIMMOMATIC_OUTPUT          } from '../../modules/local/parse_trimmomatic_output'
include { REDUCE_READS                      } from '../../modules/local/reduce_reads'

workflow TRIMMOMATIC {
    take:
    reads // channel: [ meta(id), [ reads ] ]

    main:
    ch_reports = Channel.empty()
    ch_versions = Channel.empty()

    GUNZIP (
        reads
    )

    REMOVE_POLY_GS (
        GUNZIP.out.unzipped_reads,
        params.xg
    )

    TRIMMOMATIC_MODULE (
        REMOVE_POLY_GS.out.nog_reads,
        Channel.fromPath(params.illuminaclip).first()
    )

    PARSE_TRIMMOMATIC_OUTPUT (
        TRIMMOMATIC_MODULE.out.log,
        params.min_reads
    )

    PARSE_TRIMMOMATIC_OUTPUT.out.csv
        .map {
            meta, csv ->
            [ meta ] + csv.splitCsv().collect()
        }
        .multiMap {
            meta, both_surviving, max_read_len ->
            both_surviving: [ meta, both_surviving[0].toInteger() ]
            max_read_len: [ meta, max_read_len[0].toInteger() ]
        }
        .set { ch_output }

    TRIMMOMATIC_MODULE.out.trimmed_reads
        .join(ch_output.both_surviving)
        .branch {
            meta, reads, both_surviving ->
            reduce: both_surviving > params.read_cutoff
            skip: true
        }
        .set { ch_reduce }

    REDUCE_READS (
        ch_reduce.reduce
            .map {
                meta, reads, both_surviving ->
                [ meta, reads ]
            },
        true,
        params.read_cutoff
    )

    ch_reduce.skip
        .map {
            meta, reads, both_surviving ->
            [ meta, reads ]
        }
        .set { ch_skipped_reads }

    // Collect reports
    ch_reports = ch_reports.concat(PARSE_TRIMMOMATIC_OUTPUT.out.report)

    // Collect versions
    ch_versions = ch_versions.mix(GUNZIP.out.versions)
    ch_versions = ch_versions.mix(REMOVE_POLY_GS.out.versions)
    ch_versions = ch_versions.mix(TRIMMOMATIC_MODULE.out.versions)
    ch_versions = ch_versions.mix(PARSE_TRIMMOMATIC_OUTPUT.out.versions)
    ch_versions = ch_versions.mix(REDUCE_READS.out.versions)

    emit:
    reads = REDUCE_READS.out.reduced_reads.mix(ch_skipped_reads)
    max_read_len = ch_output.max_read_len
    reports = ch_reports
    versions = ch_versions // channel: [ versions.yml ]
}
