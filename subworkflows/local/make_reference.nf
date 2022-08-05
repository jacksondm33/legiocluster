include { BWA_INDEX          } from '../../modules/local/bwa_index'
include { SAMTOOLS_FAIDX     } from '../../modules/local/samtools_faidx'
include { CREATE_SNP_CONS_FA } from '../../modules/local/create_snp_cons_fa'
include { TOUCH              } from '../../modules/local/touch'

workflow MAKE_REFERENCE {
    take:
    fasta // channel: [ meta(ref), fasta ]

    main:
    ch_reports = Channel.empty()
    ch_versions = Channel.empty()

    // Run bwa index
    BWA_INDEX (
        fasta
    )

    // Run samtools faidx
    SAMTOOLS_FAIDX (
        fasta
    )

    // Create SNP consensus (fasta)
    CREATE_SNP_CONS_FA (
        fasta
    )

    TOUCH (
        fasta
            .map {
                meta, fasta ->
                [ meta, "mutations_matrix.csv" ]
            }
    )

    emit:
    fasta = fasta
    bwa = BWA_INDEX.out.index
    fai = SAMTOOLS_FAIDX.out.fai
    snp_cons = CREATE_SNP_CONS_FA.out.snp_cons
    mutations_matrix = TOUCH.out.touch
    reports = ch_reports
    versions = ch_versions // channel: [ versions.yml ]
}
