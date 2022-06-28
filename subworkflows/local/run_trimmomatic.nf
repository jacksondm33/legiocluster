//
// Check input samplesheet and get read channels
//

include { TRIMMOMATIC } from '../../modules/nf-core/modules/trimmomatic/main'
include { REMOVE_POLY_GS } from '../../modules/local/remove_poly_gs'
include { PARSE_TRIMMOMATIC_LOG } from '../../modules/local/parse_trimmomatic_log'

workflow RUN_TRIMMOMATIC {
    take:
    reads

    main:
    ch_versions = Channel.empty()

    REMOVE_POLY_GS (
        reads
    )
    ch_versions = ch_versions.mix(REMOVE_POLY_GS.out.versions)

    TRIMMOMATIC (
        REMOVE_POLY_GS.out.reads_nog
    )
    ch_verions = ch_versions.mix(TRIMMOMATIC.out.versions)

    PARSE_TRIMMOMATIC_LOG (
        TRIMMOMATIC.out.log
    )
    ch_verions = ch_versions.mix(PARSE_TRIMMOMATIC_LOG.out.versions)

    (both_surviving, max_read_len) = PARSE_TRIMMOMATIC_LOG.out.summary.splitText().collect()
    if (both_surviving < params.min_reads) {
        error "Not enough reads surviving after Trimmomatic."
    }

    READ_REDUCER (
        TRIMMOMATIC.out.trimmed_reads
        both_surviving
    )
    ch_verions = ch_versions.mix(READ_REDUCER.out.versions)

    emit:
    remove_poly_gs_log = REMOVE_POLY_GS.out.log
    log = TRIMMOMATIC.out.log
    report = PARSE_TRIMMOMATIC_LOG.out.report
    reads = READ_REDUCER.out ? READ_REDUCER.out.reads : TRIMMOMATIC.out.trimmed_reads
    versions = ch_versions // channel: [ versions.yml ]
}
