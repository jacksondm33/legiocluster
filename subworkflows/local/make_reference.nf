include { BWA_INDEX        } from '../../modules/local/bwa_index'
include { SAMTOOLS_FAIDX   } from '../../modules/local/samtools_faidx'
include { MAKE_SNP_CONS_FA } from '../../modules/local/make_snp_cons_fa'
include { TOUCH            } from '../../modules/local/touch'

workflow MAKE_REFERENCE {
    take:
    fasta // channel: [ meta(ref), fasta ]

    main:
    ch_versions = Channel.empty()

    // Run bwa index
    BWA_INDEX (
        fasta
    )

    // Run samtools faidx
    SAMTOOLS_FAIDX (
        fasta
    )

    // Make SNP consensus (fasta)
    MAKE_SNP_CONS_FA (
        fasta
    )

    TOUCH (
        fasta
            .map {
                meta, fasta ->
                [ meta, "mutations_matrix.csv" ]
            }
    )

    // Collect versions
    ch_versions = ch_versions.mix(BWA_INDEX.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    ch_versions = ch_versions.mix(MAKE_SNP_CONS_FA.out.versions)

    emit:
    fasta = fasta
    bwa = BWA_INDEX.out.bwa
    fai = SAMTOOLS_FAIDX.out.fai
    snp_cons = MAKE_SNP_CONS_FA.out.snp_cons
    mutations_matrix = TOUCH.out.touch
    versions = ch_versions // channel: [ versions.yml ]
}
