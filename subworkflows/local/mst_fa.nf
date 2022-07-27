include { CREATE_SNP_CONS } from '../../modules/local/create_snp_cons'

workflow MST_FA {
    take:
    fasta // channel: [ val(meta), fasta ]
    vcfs

    main:
    ch_reports = Channel.empty()
    ch_versions = Channel.empty()

    CREATE_SNP_CONS (
        fasta.join(vcfs)
    )

    // Collect reports

    // Collect versions
    ch_versions = ch_versions.mix(CREATE_SNP_CONS.out.versions)

    emit:
    bases = CREATE_SNP_CONS.out.bases
    reports = ch_reports
    versions = ch_versions // channel: [ versions.yml ]
}
