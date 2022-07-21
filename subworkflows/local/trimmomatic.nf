include { GUNZIP                            } from '../../modules/local/gunzip'
include { REMOVE_POLY_GS                    } from '../../modules/local/remove_poly_gs'
include { TRIMMOMATIC as TRIMMOMATIC_MODULE } from '../../modules/local/trimmomatic'
include { PARSE_TRIMMOMATIC_LOG             } from '../../modules/local/parse_trimmomatic_log'
include { REDUCE_READS                      } from '../../modules/local/reduce_reads'

workflow TRIMMOMATIC {
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

    TRIMMOMATIC_MODULE (
        REMOVE_POLY_GS.out.nog_reads,
        Channel.fromPath(params.illuminaclip).first()
    )

    PARSE_TRIMMOMATIC_LOG (
        TRIMMOMATIC_MODULE.out.log
    )

    PARSE_TRIMMOMATIC_LOG.out.csv
        .map {
            meta, csv ->
            [ meta ] + csv.splitCsv().collect { it[0] }
        }
        .multiMap {
            meta, both_surviving, max_read_len ->
            both_surviving: [meta, both_surviving.toInteger()]
            max_read_len: [meta, max_read_len.toInteger()]
        }
        .set { ch_output }


    REDUCE_READS (
        TRIMMOMATIC_MODULE.out.trimmed_reads.join(ch_output.both_surviving)
    )

    // Collect reports
    ch_reports = ch_reports.concat(PARSE_TRIMMOMATIC_LOG.out.report)

    // Collect versions
    ch_versions = ch_versions.mix(GUNZIP.out.versions)
    ch_versions = ch_versions.mix(REMOVE_POLY_GS.out.versions)
    ch_versions = ch_versions.mix(TRIMMOMATIC_MODULE.out.versions)
    ch_versions = ch_versions.mix(PARSE_TRIMMOMATIC_LOG.out.versions)
    ch_versions = ch_versions.mix(REDUCE_READS.out.versions)

    emit:
    reads = REDUCE_READS.out.reduced_reads
    max_read_len = ch_output.max_read_len
    reports = ch_reports
    versions = ch_versions // channel: [ versions.yml ]
}
