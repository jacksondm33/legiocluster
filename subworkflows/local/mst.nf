include { CLEANUP_FREEBAYES } from '../../modules/local/cleanup_freebayes'

workflow MST {
    take:
    fasta // channel: [ val(meta), fasta ]
    mpileup
    vcf

    main:
    ch_reports = Channel.empty()
    ch_versions = Channel.empty()

    CLEANUP_FREEBAYES (
        fasta.join(mpileup).join(vcf),
        false
    )

    // COMPARE_SNP (
    // )

    // MAKE_MST (
    // )

    // Collect reports
    // ch_reports = ch_reports.concat(PARSE_VCF_OUTPUT.out.report)

    // Collect versions
    // ch_versions = ch_versions.mix(CLEANUP_FREEBAYES.out.versions)

    emit:
    reports = ch_reports
    versions = ch_versions // channel: [ versions.yml ]
}
