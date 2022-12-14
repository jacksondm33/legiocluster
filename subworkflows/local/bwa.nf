include { BWA_MEM                                 } from '../../modules/local/bwa_mem'
include { SAMTOOLS_SORT                           } from '../../modules/local/samtools_sort'
include { SAMTOOLS_INDEX                          } from '../../modules/local/samtools_index'
include { PICARD_MARKDUPLICATES                   } from '../../modules/local/picard_markduplicates'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MARKED } from '../../modules/local/samtools_index'
include { BCFTOOLS_MPILEUP                        } from '../../modules/local/bcftools_mpileup'
include { BCFTOOLS_VIEW                           } from '../../modules/local/bcftools_view'
include { SAMTOOLS_FLAGSTAT                       } from '../../modules/local/samtools_flagstat'
include { SAMTOOLS_IDXSTATS                       } from '../../modules/local/samtools_idxstats'
include { SAMTOOLS_DEPTH                          } from '../../modules/local/samtools_depth'
include { PARSE_BWA_OUTPUT                        } from '../../modules/local/parse_bwa_output'

workflow BWA {
    take:
    reads            // channel: [ meta(id, ref), [ reads ]        ]
    fasta            // channel: [ meta(id, ref), fasta            ]
    bwa              // channel: [ meta(id, ref), bwa              ]
    fai              // channel: [ meta(id, ref), fai              ]
    mapped_threshold // channel: [ meta(id, ref), mapped_threshold ]

    main:
    ch_reports = Channel.empty()
    ch_versions = Channel.empty()

    BWA_MEM (
        reads.join(bwa)
    )

    SAMTOOLS_SORT (
        BWA_MEM.out.sam
    )

    SAMTOOLS_INDEX (
        SAMTOOLS_SORT.out.bam
    )

    PICARD_MARKDUPLICATES (
        SAMTOOLS_SORT.out.bam
    )

    SAMTOOLS_INDEX_MARKED (
        PICARD_MARKDUPLICATES.out.marked_bam
    )

    BCFTOOLS_MPILEUP (
        PICARD_MARKDUPLICATES.out.marked_bam
            .join(fasta)
            .join(fai)
    )

    BCFTOOLS_VIEW (
        BCFTOOLS_MPILEUP.out.bcf
    )

    SAMTOOLS_FLAGSTAT (
        PICARD_MARKDUPLICATES.out.marked_bam.join(SAMTOOLS_INDEX_MARKED.out.bai)
    )

    SAMTOOLS_IDXSTATS (
        PICARD_MARKDUPLICATES.out.marked_bam.join(SAMTOOLS_INDEX_MARKED.out.bai)
    )

    SAMTOOLS_DEPTH (
        PICARD_MARKDUPLICATES.out.marked_bam
    )

    PARSE_BWA_OUTPUT (
        fasta
            .join(BWA_MEM.out.sam)
            .join(SAMTOOLS_FLAGSTAT.out.flagstat)
            .join(SAMTOOLS_IDXSTATS.out.idxstats)
            .join(mapped_threshold)
    )

    PARSE_BWA_OUTPUT.out.csv
        .map {
            meta, csv ->
            [ meta ] + csv.splitCsv().collect()
        }
        .map {
            meta, percent_mapped ->
            [ meta, percent_mapped[0].toFloat() ]
        }
        .set { ch_percent_mapped }

    // Collect reports
    ch_reports = ch_reports.concat(PARSE_BWA_OUTPUT.out.report)

    // Collect versions
    ch_versions = ch_versions.mix(BWA_MEM.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)
    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_MARKED.out.versions)
    ch_versions = ch_versions.mix(BCFTOOLS_MPILEUP.out.versions)
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_DEPTH.out.versions)
    ch_versions = ch_versions.mix(PARSE_BWA_OUTPUT.out.versions)

    emit:
    percent_mapped = ch_percent_mapped
    depth = SAMTOOLS_DEPTH.out.depth
    bam = PICARD_MARKDUPLICATES.out.marked_bam
    mpileup = BCFTOOLS_VIEW.out.vcf
    reports = ch_reports
    versions = ch_versions // channel: [ versions.yml ]
}
