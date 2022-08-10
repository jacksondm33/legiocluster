include { KRAKEN as KRAKEN_MODULE } from '../../modules/local/kraken'
include { PARSE_KRAKEN_OUTPUT     } from '../../modules/local/parse_kraken_output'

workflow KRAKEN {
    take:
    fasta // channel: [ meta(id, ref), fasta ]

    main:
    ch_reports = Channel.empty()
    ch_versions = Channel.empty()

    KRAKEN_MODULE (
        fasta,
        Channel.fromPath(params.kraken_db, type: 'dir').first()
    )

    PARSE_KRAKEN_OUTPUT (
        KRAKEN_MODULE.out.results.join(fasta),
        params.genomes.get(params.genome).binomial.tokenize()[0]
    )

    // Collect reports
    ch_reports = ch_reports.concat(PARSE_KRAKEN_OUTPUT.out.report)

    // Collect versions
    ch_versions = ch_versions.mix(KRAKEN_MODULE.out.versions)
    ch_versions = ch_versions.mix(PARSE_KRAKEN_OUTPUT.out.versions)

    emit:
    good_contigs = PARSE_KRAKEN_OUTPUT.out.good_contigs
    bad_contigs = PARSE_KRAKEN_OUTPUT.out.bad_contigs
    reports = ch_reports
    versions = ch_versions // channel: [ versions.yml ]
}
