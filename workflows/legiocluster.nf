/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowLegiocluster.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.references, params.genomes_mash, params.illuminaclip, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Create input and references channels
ch_input = file(params.input)
ch_references = file(params.references)
ch_genomes_mash = file(params.genomes_mash)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Modules
include { MASH_SKETCH                 } from '../modules/local/mash_sketch'
include { CREATE_SNP_CONS_FA          } from '../modules/local/create_snp_cons_fa'
include { CREATE_SNP_CONS             } from '../modules/local/create_snp_cons'
include { COMPARE_SNPS                } from '../modules/local/compare_snps'
include { CREATE_REPORT               } from '../modules/local/create_report'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/local/custom_dumpsoftwareversions'
include { MULTIQC                     } from '../modules/local/multiqc'

// Subworkflows
include { CHECK_INPUT } from '../subworkflows/local/check_input'
include { TRIMMOMATIC } from '../subworkflows/local/trimmomatic'
include { FASTQC      } from '../subworkflows/local/fastqc'
include { MASH_FQ     } from '../subworkflows/local/mash_fq'
include { SPADES      } from '../subworkflows/local/spades'
include { MASH_FA     } from '../subworkflows/local/mash_fa'
include { BWA_FA      } from '../subworkflows/local/bwa_fa'
include { BWA         } from '../subworkflows/local/bwa'
include { QUAST       } from '../subworkflows/local/quast'
include { QUALIMAP    } from '../subworkflows/local/qualimap'
include { FREEBAYES   } from '../subworkflows/local/freebayes'
include { MAKE_MST    } from '../subworkflows/local/make_mst'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow LEGIOCLUSTER {

    ch_reports = Channel.empty()
    ch_versions = Channel.empty()

    // Check input
    CHECK_INPUT (
        ch_input,
        ch_references
    )

    // Run trimmomatic
    TRIMMOMATIC (
        CHECK_INPUT.out.reads
    )

    // Run fastqc
    FASTQC (
        CHECK_INPUT.out.reads,
        TRIMMOMATIC.out.reads
    )

    // Run mash fq
    // Combine trimmomatic reads
    // with genomes mash file
    TRIMMOMATIC.out.reads
        .combine(ch_genomes_mash)
        .multiMap {
            meta, reads, mash ->
            reads: [ meta, reads ]
            mash:  [ meta, mash  ]
        }
        .set { ch_mash_fq_in }
    MASH_FQ (
        ch_mash_fq_in.reads,
        ch_mash_fq_in.mash
    )

    // Run spades
    SPADES (
        TRIMMOMATIC.out.reads,
        TRIMMOMATIC.out.max_read_len
    )

    // Run mash fa
    // Combine spades filtered contigs
    // with list of reference fastas
    SPADES.out.filtered_contigs
        .combine(
            CHECK_INPUT.out.fasta
                .collect { it[1] }
        )
        .multiMap {
            meta, fasta, fastas ->
            fasta:  [ meta, fasta  ]
            fastas: [ meta, fastas ]
        }
        .set { ch_mash_fa_in }
    MASH_FA (
        ch_mash_fa_in.fasta,
        ch_mash_fa_in.fastas
    )

    // Create bwa reads, fasta channel [ meta(id, ref), reads, fasta ]
    TRIMMOMATIC.out.reads
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
        .cross(CHECK_INPUT.out.fasta) { it[0].ref }
        .map { it[0] + it[1][1..-1] }
        .set { ch_bwa_in }

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

    // CREATE_REPORT (
    //     ch_reports
    //         .map {
    //             meta, report ->
    //             [ meta - [ref: meta.ref], report ]
    //         }
    //         .groupTuple()
    //         .join(CHECK_INPUT.out.reads),
    //     params.genome
    // )

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    // MultiQC
    workflow_summary    = WorkflowLegiocluster.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.results.collect { it[1] }.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
