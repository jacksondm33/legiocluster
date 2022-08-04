// Modules
include { MASH_SKETCH                 } from '../modules/local/mash_sketch'
include { CREATE_SNP_CONS_FA          } from '../modules/local/create_snp_cons_fa'
include { CREATE_SNP_CONS             } from '../modules/local/create_snp_cons'
include { COMPARE_SNPS                } from '../modules/local/compare_snps'
include { CHECK_REF_QUAL              } from '../modules/local/check_ref_qual'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/local/custom_dumpsoftwareversions'
include { MULTIQC                     } from '../modules/local/multiqc'

// Subworkflows
include { MASH_FA     } from '../subworkflows/local/mash_fa'
include { BWA_FA      } from '../subworkflows/local/bwa_fa'
include { BWA         } from '../subworkflows/local/bwa'
include { QUAST       } from '../subworkflows/local/quast'
include { QUALIMAP    } from '../subworkflows/local/qualimap'
include { FREEBAYES   } from '../subworkflows/local/freebayes'
include { MAKE_MST    } from '../subworkflows/local/make_mst'

workflow LEGIOCLUSTER_MAIN {
    take:
    ch_reads            // channel: [ meta(id), [ reads ] ]
    ch_contigs          // channel: [ meta(id), contigs ]
    ch_filtered_contigs // channel: [ meta(id), filtered_contigs ]
    ch_fastas           // channel: [ [ meta(ref), fasta ] ]
    ch_bwa_index        // channel: [ [ meta(ref), bwa_index ] ]
    ch_bwa_fai          // channel: [ [ meta(ref), bwa_fai ] ]
    ch_reports          // channel: [ report.txt ]
    ch_versions         // channel: [ versions.yml ]

    main:

    ch_fastas
        .flatMap()
        .set { ch_fasta }

    ch_filtered_contigs.dump(tag: 'contigs')
    ch_fastas.dump(tag: 'fastas')
    ch_fasta.dump(tag: 'fasta')

    // Mash sketch
    MASH_SKETCH (
        ch_fasta
            .collect { it[1] }
            .map { [ [:], it ] },
        false,
        false
    )

    // Run mash fa
    MASH_FA (
        ch_filtered_contigs,
        MASH_SKETCH.out.mash.map { it[1] }
    )

    // Create bwa reads, fasta channel [ meta(id, ref), reads, fasta ]
    ch_reads
        .combine(
            MASH_FA.out.fastas
                .map {
                    meta, fastas ->
                    [ meta, meta.set_ref != '' ? [ meta.set_ref ] : fastas ]
                }
                .transpose(),
            by: 0)
        .map {
            meta, reads, fasta ->
            [ meta + [ref: fasta], reads ]
        }
        .cross(ch_fasta) { it[0].ref }
        .map { it[0] + [ it[1][1] ] }
        .set { ch_bwa_reads_fasta }

    ch_bwa_reads_fasta
        .map {
            meta, reads, fasta ->
            [ [ref: meta.ref], fasta ]
        }

    // Run bwa fa
    BWA_FA (
        ch_bwa_reads_fasta
            .map {
                meta, reads, fasta ->
                [ [ref: meta.ref], fasta ]
            }
            .unique()
    )

    // Create bwa input channels
    ch_bwa_reads_fasta
        .cross(
            BWA_FA.out.index
                .join(BWA_FA.out.fai)
        ) { it[0].ref }
        .map { it[0] + [ it[1][1], it[1][2] ] }
        .multiMap {
            meta, reads, fasta, index, fai ->
            reads: [ meta, reads ]
            fasta: [ meta, fasta ]
            index: [ meta, index ]
            fai: [ meta, fai ]
            mapped_threshold: [ meta, meta.set_ref != 'NO_FILE' ? 0 : meta.make_ref == 'true' ? 100 : params.mapped_threshold ]
        }
        .set { ch_bwa_input }

    // Run bwa
    BWA (
        ch_bwa_input.reads,
        ch_bwa_input.fasta,
        ch_bwa_input.index,
        ch_bwa_input.fai,
        ch_bwa_input.mapped_threshold
    )

    BWA.out.percent_mapped
        .join(BWA.out.depth)
        .join(BWA.out.bam)
        .join(BWA.out.mpileup)
        .join(ch_bwa_input.fasta)
        .join(ch_bwa_input.fai)
        .join(ch_bwa_input.mapped_threshold)
        .map {
            meta, percent_mapped, depth, bam, mpileup, fasta, fai, mapped_threshold ->
            [ meta - [ref: meta.ref], percent_mapped, depth, bam, mpileup, fasta, fai, mapped_threshold ]
        }
        .join(ch_contigs)
        .map {
            meta, percent_mapped, depth, bam, mpileup, fasta, fai, mapped_threshold, contigs ->
            [ meta, [ percent_mapped, depth, bam, mpileup, fasta, fai, mapped_threshold, contigs ] ]
        }
        .groupTuple()
        .map {
            meta, output ->
            [ meta ] + output.max { it[0] }
        }
        .map {
            meta, percent_mapped, depth, bam, mpileup, fasta, fai, mapped_threshold, contigs ->
            [ meta + [ref: fasta.baseName], percent_mapped, depth, bam, mpileup, fasta, fai, mapped_threshold, contigs ]
        }
        .multiMap {
            meta, percent_mapped, depth, bam, mpileup, fasta, fai, mapped_threshold, contigs ->
            percent_mapped: [ meta, percent_mapped ]
            depth: [ meta, depth ]
            bam: [ meta, bam ]
            mpileup: [ meta, mpileup ]
            fasta: [ meta, fasta ]
            fai: [ meta, fai ]
            contigs: [ meta, contigs ]
            max_no_ns:          [ meta, meta.set_ref != 'NO_FILE' ? 999999999 : meta.make_ref == 'true' ? 999999999 : params.max_no_ns          ]
            max_no_gaps:        [ meta, meta.set_ref != 'NO_FILE' ? 999999999 : meta.make_ref == 'true' ? 999999999 : params.max_no_gaps        ]
            snp_threshold:      [ meta, meta.set_ref != 'NO_FILE' ? 999999999 : meta.make_ref == 'true' ? 0         : params.snp_threshold      ]
            min_percent_mapped: [ meta, meta.set_ref != 'NO_FILE' ? 0         : meta.make_ref == 'true' ? 0         : params.min_percent_mapped ]
            mapped_threshold:   [ meta, mapped_threshold ]
        }
        .set { ch_bwa_output }

    // Run quast
    QUAST (
        ch_bwa_output.contigs,
        ch_bwa_output.fasta,
        ch_bwa_output.depth,
        ch_bwa_output.percent_mapped,
        ch_bwa_output.max_no_ns,
        ch_bwa_output.max_no_gaps,
        ch_bwa_output.min_percent_mapped,
        ch_bwa_output.mapped_threshold
    )

    // Run qualimap
    QUALIMAP (
        ch_bwa_output.bam
    )

    // Run freebayes
    FREEBAYES (
        ch_bwa_output.bam,
        ch_bwa_output.fasta,
        ch_bwa_output.fai,
        ch_bwa_output.snp_threshold,
        QUAST.out.depth_mean,
        QUAST.out.depth_sd
    )

    FREEBAYES.out.mutations
        .join(ch_bwa_output.percent_mapped)
        .join(ch_bwa_output.snp_threshold)
        .join(ch_bwa_output.mapped_threshold)
        .branch {
            meta, mutations, percent_mapped, snp_threshold, mapped_threshold ->
            close: (mutations[0] < snp_threshold) && (percent_mapped > mapped_threshold)
            distant: true
        }
        .set { ch_distances }

    // Freebayes close distances channel
    ch_distances.close
        .join(ch_bwa_output.fasta)
        .join(ch_bwa_output.mpileup)
        .join(FREEBAYES.out.vcf)
        .multiMap {
            meta, mutations, percent_mapped, snp_threshold, mapped_threshold, fasta, mpileup, freebayes ->
            fasta:     [ meta, fasta ]
            mpileup:   [ meta, mpileup ]
            freebayes: [ meta, freebayes ]
        }
        .set { ch_freebayes_close }

    ch_distances.distant
        .join(ch_bwa_output.contigs)
        .multiMap {
            meta, mutations, percent_mapped, snp_threshold, mapped_threshold, contigs ->
            contigs: [ meta, contigs ]
        }
        .set { ch_freebayes_distant }

    CHECK_REF_QUAL (
        ch_freebayes_distant.contigs,
        params.med_genome_len
    )

    CHECK_REF_QUAL.out.csv
        .map {
            meta, csv ->
            [ meta ] + csv.splitCsv().collect()
        }
        .map {
            meta, passed_qc ->
            [ meta, passed_qc[0].toBoolean() ]
        }
        .branch {
            meta, passed_qc ->
            make_ref: passed_qc
        }
        .set { ch_ref }

    // Collect reports
    ch_reports = ch_reports.concat(MASH_FA.out.reports)
    ch_reports = ch_reports.concat(BWA.out.reports)
    ch_reports = ch_reports.concat(QUAST.out.reports)
    ch_reports = ch_reports.concat(QUALIMAP.out.reports)
    ch_reports = ch_reports.concat(FREEBAYES.out.reports)

    // Collect versions
    ch_versions = ch_versions.mix(MASH_FA.out.versions)
    ch_versions = ch_versions.mix(BWA_FA.out.versions)
    ch_versions = ch_versions.mix(BWA.out.versions)
    ch_versions = ch_versions.mix(QUAST.out.versions)
    ch_versions = ch_versions.mix(QUALIMAP.out.versions)
    ch_versions = ch_versions.mix(FREEBAYES.out.versions)

    emit:
    Channel.empty()     // channel: [ meta(id), [ reads ] ]
    Channel.empty()     // channel: [ meta(id), contigs ]
    Channel.empty()     // channel: [ meta(id), filtered_contigs ]
    Channel.empty()     // channel: [ meta(ref), fasta ]
    BWA_FA.out.index    // channel: [ meta(ref), bwa_index ]
    BWA_FA.out.fai      // channel: [ meta(ref), bwa_fai ]
    ch_reports          // channel: [ report.txt ]
    ch_versions         // channel: [ versions.yml ]
}
