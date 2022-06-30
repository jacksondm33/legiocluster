include { GUNZIP } from '../../modules/nf-core/modules/gunzip/main'
include { REMOVE_POLY_GS } from '../../modules/local/remove_poly_gs'
include { TRIMMOMATIC } from '../../modules/local/trimmomatic'
include { PARSE_TRIMMOMATIC_LOG } from '../../modules/local/parse_trimmomatic_log'
include { READ_REDUCER } from '../../modules/local/read_reducer'

workflow RUN_TRIMMOMATIC {
    take:
    reads

    main:
    ch_versions = Channel.empty()

    GUNZIP (
        reads.transpose()
    )

    REMOVE_POLY_GS (
        GUNZIP.out.gunzip.groupTuple()
    )

    TRIMMOMATIC (
        REMOVE_POLY_GS.out.nog_reads,
        Channel.fromPath(params.nextera_pe)
    )

    PARSE_TRIMMOMATIC_LOG (
        TRIMMOMATIC.out.log
    )

    READ_REDUCER (
        TRIMMOMATIC.out.trimmed_reads,
        PARSE_TRIMMOMATIC_LOG.out.both_surviving.map { meta, both_surviving -> both_surviving.toInteger() }
    )

    // Collect versions
    ch_versions = ch_versions.mix(GUNZIP.out.versions)
    ch_versions = ch_versions.mix(REMOVE_POLY_GS.out.versions)
    ch_versions = ch_versions.mix(TRIMMOMATIC.out.versions)
    ch_versions = ch_versions.mix(PARSE_TRIMMOMATIC_LOG.out.versions)
    ch_versions = ch_versions.mix(READ_REDUCER.out.versions)

    emit:
    reads = READ_REDUCER.out.reduced_reads
    max_read_len = PARSE_TRIMMOMATIC_LOG.out.max_read_len
    versions = ch_versions // channel: [ versions.yml ]
}
