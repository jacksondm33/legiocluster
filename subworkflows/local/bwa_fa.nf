include { BWA_INDEX      } from '../../modules/local/bwa_index'
include { SAMTOOLS_FAIDX } from '../../modules/local/samtools_faidx'

workflow BWA_FA {
    take:
    fasta // channel: [ meta(ref), fasta ]

    main:
    ch_versions = Channel.empty()

    BWA_INDEX (
        fasta
    )

    SAMTOOLS_FAIDX (
        fasta
    )

    // Collect versions
    ch_versions = ch_versions.mix(BWA_INDEX.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    emit:
    bwa = BWA_INDEX.out.index
    fai = SAMTOOLS_FAIDX.out.fai
    versions = ch_versions // channel: [ versions.yml ]
}
