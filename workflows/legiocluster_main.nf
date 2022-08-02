// // Modules
include { MASH_SKETCH } from '../modules/local/mash_sketch'
// include { MASH_SKETCH as MASH_SKETCH_STRAINS } from '../modules/local/mash_sketch'
// include { CREATE_SNP_CONS_FA                 } from '../modules/local/create_snp_cons_fa'
// include { CREATE_SNP_CONS                    } from '../modules/local/create_snp_cons'
// include { COMPARE_SNPS                       } from '../modules/local/compare_snps'
// include { CREATE_REPORT                      } from '../modules/local/create_report'
// include { CUSTOM_DUMPSOFTWAREVERSIONS        } from '../modules/local/custom_dumpsoftwareversions'
// include { MULTIQC                            } from '../modules/local/multiqc'

// // Subworkflows
// include { CHECK_INPUT } from '../subworkflows/local/check_input'
// include { TRIMMOMATIC } from '../subworkflows/local/trimmomatic'
// include { FASTQC      } from '../subworkflows/local/fastqc'
// include { MASH_FQ     } from '../subworkflows/local/mash_fq'
// include { SPADES      } from '../subworkflows/local/spades'
// include { MASH_FA     } from '../subworkflows/local/mash_fa'
// include { BWA_FA      } from '../subworkflows/local/bwa_fa'
// include { BWA         } from '../subworkflows/local/bwa'
// include { QUAST       } from '../subworkflows/local/quast'
// include { QUALIMAP    } from '../subworkflows/local/qualimap'
// include { FREEBAYES   } from '../subworkflows/local/freebayes'
// include { MAKE_MST    } from '../subworkflows/local/make_mst'

workflow LEGIOCLUSTER_MAIN {
    take:
    reads
    contigs
    filtered_contigs
    fastas
    reports
    versions

    main:
    // reads.dump(tag: 'reads')
    contigs.dump(tag: 'contigs')
    fastas.dump(tag: 'fastas')

    // Mash sketch
    MASH_SKETCH (
        fastas,
        false,
        false
    )

    // // Run mash fa
    // MASH_FA (
    //     SPADES.out.filtered_contigs,
    //     MASH_SKETCH_STRAINS.out.mash // TODO: This has an id now
    // )

    // // Create bwa reads, fasta channel [ meta(id, ref), reads, fasta ]
    // TRIMMOMATIC.out.reads
    //     .combine(
    //         MASH_FA.out.fastas
    //             .map {
    //                 meta, fastas ->
    //                 [ meta, meta.set_ref != 'NO_FILE' ? [ meta.set_ref ] : fastas ]
    //             }
    //             .transpose(),
    //         by: 0)
    //     .map {
    //         meta, reads, fasta ->
    //         [ meta + [ref: file(fasta).baseName], reads ]
    //     }
    //     .cross(ch_strain_fasta) { it[0].ref } // TODO: Use spades filtered_contigs
    //     .map { it[0] + [ it[1][1] ] }
    //     .set { ch_bwa_reads_fasta }

    // // Run bwa fa
    // BWA_FA (
    //     ch_bwa_reads_fasta
    //         .map {
    //             meta, reads, fasta ->
    //             [ [ref: meta.ref], fasta ]
    //         }
    //         .unique()
    // )

    // // Create bwa input channels
    // ch_bwa_reads_fasta
    //     .cross(
    //         BWA_FA.out.index
    //             .join(BWA_FA.out.fai)
    //     ) { it[0].ref }
    //     .map { it[0] + [ it[1][1], it[1][2] ] }
    //     .multiMap {
    //         meta, reads, fasta, index, fai ->
    //         reads: [ meta, reads ]
    //         fasta: [ meta, fasta ]
    //         index: [ meta, index ]
    //         fai: [ meta, fai ]
    //         mapped_threshold: [ meta, meta.set_ref != 'NO_FILE' ? 0 : meta.make_ref == 'true' ? 100 : params.mapped_threshold ]
    //     }
    //     .set { ch_bwa_input }

    // // Run bwa
    // BWA (
    //     ch_bwa_input.reads,
    //     ch_bwa_input.fasta,
    //     ch_bwa_input.index,
    //     ch_bwa_input.fai,
    //     ch_bwa_input.mapped_threshold
    // )

    // // TODO: Make sure selected reference is actually valid

    // BWA.out.percent_mapped
    //     .join(BWA.out.depth)
    //     .join(BWA.out.bam)
    //     .join(BWA.out.mpileup)
    //     .join(ch_bwa_input.fasta)
    //     .join(ch_bwa_input.fai)
    //     .join(ch_bwa_input.mapped_threshold)
    //     .map {
    //         meta, percent_mapped, depth, bam, mpileup, fasta, fai, mapped_threshold ->
    //         [ meta - [ref: meta.ref], percent_mapped, depth, bam, mpileup, fasta, fai, mapped_threshold ]
    //     }
    //     .join(SPADES.out.contigs)
    //     .map {
    //         meta, percent_mapped, depth, bam, mpileup, fasta, fai, mapped_threshold, contigs ->
    //         [ meta, [ percent_mapped, depth, bam, mpileup, fasta, fai, mapped_threshold, contigs ] ]
    //     }
    //     .groupTuple()
    //     .map {
    //         meta, output ->
    //         [ meta ] + output.max { it[0] }
    //     }
    //     .map {
    //         meta, percent_mapped, depth, bam, mpileup, fasta, fai, mapped_threshold, contigs ->
    //         [ meta + [ref: fasta.baseName], percent_mapped, depth, bam, mpileup, fasta, fai, mapped_threshold, contigs ]
    //     }
    //     .multiMap {
    //         meta, percent_mapped, depth, bam, mpileup, fasta, fai, mapped_threshold, contigs ->
    //         percent_mapped: [ meta, percent_mapped ]
    //         depth: [ meta, depth ]
    //         bam: [ meta, bam ]
    //         mpileup: [ meta, mpileup ]
    //         fasta: [ meta, fasta ]
    //         fai: [ meta, fai ]
    //         contigs: [ meta, contigs ]
    //         max_no_ns:          [ meta, meta.set_ref != 'NO_FILE' ? 999999999 : meta.make_ref == 'true' ? 999999999 : params.max_no_ns          ]
    //         max_no_gaps:        [ meta, meta.set_ref != 'NO_FILE' ? 999999999 : meta.make_ref == 'true' ? 999999999 : params.max_no_gaps        ]
    //         snp_threshold:      [ meta, meta.set_ref != 'NO_FILE' ? 999999999 : meta.make_ref == 'true' ? 0         : params.snp_threshold      ]
    //         min_percent_mapped: [ meta, meta.set_ref != 'NO_FILE' ? 0         : meta.make_ref == 'true' ? 0         : params.min_percent_mapped ]
    //         mapped_threshold:   [ meta, mapped_threshold ]
    //     }
    //     .set { ch_bwa_output }

    // // Run quast
    // QUAST (
    //     ch_bwa_output.contigs,
    //     ch_bwa_output.fasta,
    //     ch_bwa_output.depth,
    //     ch_bwa_output.percent_mapped,
    //     ch_bwa_output.max_no_ns,
    //     ch_bwa_output.max_no_gaps,
    //     ch_bwa_output.min_percent_mapped,
    //     ch_bwa_output.mapped_threshold
    // )

    // // Run qualimap
    // QUALIMAP (
    //     ch_bwa_output.bam
    // )

    // // Run freebayes
    // FREEBAYES (
    //     ch_bwa_output.bam,
    //     ch_bwa_output.fasta,
    //     ch_bwa_output.fai,
    //     ch_bwa_output.snp_threshold,
    //     QUAST.out.depth_mean,
    //     QUAST.out.depth_sd
    // )

    // FREEBAYES.out.mutations
    //     .join(ch_bwa_output.percent_mapped)
    //     .join(ch_bwa_output.snp_threshold)
    //     .join(ch_bwa_output.mapped_threshold)
    //     .branch {
    //         meta, mutations, percent_mapped, snp_threshold, mapped_threshold ->
    //         close: (mutations[0] < snp_threshold) && (percent_mapped > mapped_threshold)
    //         distant: true
    //     }
    //     .set { ch_distances }

    // // Freebayes close distances channel
    // ch_distances.close
    //     .join(ch_bwa_output.fasta)
    //     .join(ch_bwa_output.mpileup)
    //     .join(FREEBAYES.out.vcf)
    //     .multiMap {
    //         meta, mutations, percent_mapped, snp_threshold, mapped_threshold, fasta, mpileup, freebayes ->
    //         fasta:     [ meta, fasta ]
    //         mpileup:   [ meta, mpileup ]
    //         freebayes: [ meta, freebayes ]
    //     }
    //     .set { ch_freebayes_close }

    // // Create SNP consensus (fasta)
    // CREATE_SNP_CONS_FA (
    //     ch_freebayes_close.fasta
    //         .map {
    //             meta, fasta ->
    //             [ [ref: meta.ref], fasta ]
    //         }
    //         .unique()
    // )

    // // Create SNP consensus input channel
    // ch_freebayes_close.mpileup
    //     .join(ch_freebayes_close.freebayes)
    //     .cross(
    //         CREATE_SNP_CONS_FA.out.bases
    //     ) { it[0].ref }
    //     .map { it[0] + it[1][1..-1] }
    //     .set { ch_create_snp_cons_in }

    // // Create SNP consensus
    // CREATE_SNP_CONS (
    //     ch_create_snp_cons_in,
    //     false
    // )

    // // Compare SNPs input channel
    // CREATE_SNP_CONS_FA.out.snp_cons
    //     .cross(
    //         CREATE_SNP_CONS.out.snp_cons
    //     ) { it[0].ref }
    //     .groupTuple()
    //     .map { [ it[1], it[1].collect { it[1] } + [ it[0][1] ] ] }
    //     .transpose(by: 0)
    //     .map { it[0] + [ it[1] - [ it[0][1] ] ] }
    //     .set { ch_compare_snps_in }

    // COMPARE_SNPS (
    //     ch_compare_snps_in
    // )

    // COMPARE_SNPS.out.pairwise_diffs
    //     .map {
    //         meta, pairwise_diffs ->
    //         [ [ref: meta.ref], pairwise_diffs ]
    //     }
    //     .groupTuple()
    //     .set { ch_make_mst_in }

    // MAKE_MST (
    //     ch_make_mst_in
    // )

    // // Distant reference
    // // MAKE_REF (
    // //     ch_distances.distant
    // // )

    // // Collect reports
    // ch_reports = ch_reports.concat(TRIMMOMATIC.out.reports)
    // ch_reports = ch_reports.concat(FASTQC.out.reports)
    // ch_reports = ch_reports.concat(MASH_FQ.out.reports)
    // ch_reports = ch_reports.concat(SPADES.out.reports)
    // ch_reports = ch_reports.concat(MASH_FA.out.reports)
    // ch_reports = ch_reports.concat(BWA.out.reports)
    // ch_reports = ch_reports.concat(QUAST.out.reports)
    // ch_reports = ch_reports.concat(QUALIMAP.out.reports)
    // ch_reports = ch_reports.concat(FREEBAYES.out.reports)

    // CREATE_REPORT (
    //     ch_reports
    //         .map {
    //             meta, report ->
    //             [ meta - [ref: meta.ref], report ]
    //         }
    //         .groupTuple()
    //         .join(CHECK_INPUT.out.reads),
    //     params.sp_abbr
    // )

    // // Collect versions
    // ch_versions = ch_versions.mix(CHECK_INPUT.out.versions)
    // ch_versions = ch_versions.mix(TRIMMOMATIC.out.versions)
    // ch_versions = ch_versions.mix(FASTQC.out.versions)
    // ch_versions = ch_versions.mix(MASH_FQ.out.versions)
    // ch_versions = ch_versions.mix(SPADES.out.versions)
    // ch_versions = ch_versions.mix(MASH_FA.out.versions)
    // ch_versions = ch_versions.mix(BWA_FA.out.versions)
    // ch_versions = ch_versions.mix(BWA.out.versions)
    // ch_versions = ch_versions.mix(QUAST.out.versions)
    // ch_versions = ch_versions.mix(QUALIMAP.out.versions)
    // ch_versions = ch_versions.mix(FREEBAYES.out.versions)

    emit:
    Channel.empty()
    Channel.empty()
    Channel.empty()
    fastas
    reports
    versions
}
