include { CLEANUP_FREEBAYES } from '../../modules/local/cleanup_freebayes'
include { COMPARE_SNPS } from '../../modules/local/compare_snps'

workflow MST {
    take:
    mpileup
    vcf
    vcfs
    bases

    main:
    ch_reports = Channel.empty()
    ch_versions = Channel.empty()

    CLEANUP_FREEBAYES (
        mpileup.join(vcf).join(vcfs).join(bases),
        false
    )

    COMPARE_SNPS (
        vcfs.join(CLEANUP_FREEBAYES.out.csv)
    )

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
