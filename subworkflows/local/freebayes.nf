include { FREEBAYES as FREEBAYES_MODULE } from '../../modules/local/freebayes'
include { VCFFILTER                     } from '../../modules/local/vcffilter'
include { PARSE_VCF_OUTPUT              } from '../../modules/local/parse_vcf_output'

workflow FREEBAYES {
    take:
    bam   // channel: [ val(meta), bam ]
    fasta // channel: [ val(meta), fasta ]
    fai   // channel: [ val(meta), fai ]
    snp_threshold
    depth_mean
    depth_sd

    main:
    ch_reports = Channel.empty()
    ch_versions = Channel.empty()

    FREEBAYES_MODULE (
        bam.join(fasta).join(fai)
    )

    VCFFILTER (
        FREEBAYES_MODULE.out.vcf
            .join(
                depth_mean
                    .join(depth_sd)
                    .map {
                        meta, depth_mean, depth_sd ->
                        [ meta, depth_mean + 3 * depth_sd ]
                    }),
        params.qual_threshold,
        params.dp_min,
        params.qa_threshold,
        params.ao_dp_ratio
    )

    PARSE_VCF_OUTPUT (
        VCFFILTER.out.vcf.join(fasta).join(snp_threshold)
    )

    PARSE_VCF_OUTPUT.out.csv
        .map {
            meta, csv ->
            [ meta ] + csv.splitCsv().collect()
        }
        .map {
            meta, mutations ->
            [ meta, mutations.collect { it.toInteger() } ]
        }
        .set { ch_mutations }

    // Collect reports
    ch_reports = ch_reports.concat(PARSE_VCF_OUTPUT.out.report)

    // Collect versions
    ch_versions = ch_versions.mix(FREEBAYES_MODULE.out.versions)
    ch_versions = ch_versions.mix(VCFFILTER.out.versions)
    ch_versions = ch_versions.mix(PARSE_VCF_OUTPUT.out.versions)

    emit:
    mutations = ch_mutations
    vcf = VCFFILTER.out.vcf
    reports = ch_reports
    versions = ch_versions // channel: [ versions.yml ]
}
