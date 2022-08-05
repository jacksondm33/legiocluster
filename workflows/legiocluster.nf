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
ch_input = Channel.fromPath(params.input)
ch_references = Channel.fromPath(params.references)
ch_genomes_mash = Channel.fromPath(params.genomes_mash)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
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
include { CHECK_REF_QUAL              } from '../modules/local/check_ref_qual'
include { MAKE_REFERENCES             } from '../modules/local/make_references'
include { CREATE_REPORT               } from '../modules/local/create_report'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/local/custom_dumpsoftwareversions'
include { MULTIQC                     } from '../modules/local/multiqc'

// Subworkflows
include { CHECK_INPUT    } from '../subworkflows/local/check_input'
include { TRIMMOMATIC    } from '../subworkflows/local/trimmomatic'
include { FASTQC         } from '../subworkflows/local/fastqc'
include { MASH_FQ        } from '../subworkflows/local/mash_fq'
include { SPADES         } from '../subworkflows/local/spades'
include { MASH_FA        } from '../subworkflows/local/mash_fa'
include { BWA            } from '../subworkflows/local/bwa'
include { QUAST          } from '../subworkflows/local/quast'
include { QUALIMAP       } from '../subworkflows/local/qualimap'
include { FREEBAYES      } from '../subworkflows/local/freebayes'
include { MAKE_MST       } from '../subworkflows/local/make_mst'
include { MAKE_REFERENCE } from '../subworkflows/local/make_reference'

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

    // Run Trimmomatic
    TRIMMOMATIC (
        CHECK_INPUT.out.reads
    )

    // Run FastQC
    FASTQC (
        CHECK_INPUT.out.reads,
        TRIMMOMATIC.out.reads
    )

    // Mash FQ channel
    // Contains the reads for each sample
    // combined with genomes mash file
    TRIMMOMATIC.out.reads
        .combine(ch_genomes_mash)
        .multiMap {
            meta, reads, mash ->
            reads: [ meta, reads ]
            mash:  [ meta, mash  ]
        }
        .set { ch_mash_fq }

    // Run Mash FQ
    MASH_FQ (
        ch_mash_fq.reads,
        ch_mash_fq.mash
    )

    // Run SPAdes
    SPADES (
        TRIMMOMATIC.out.reads,
        TRIMMOMATIC.out.max_read_len
    )

    // Mash FA channel
    // Contains the filtered contigs for each sample
    // combined with the list of reference fastas
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
        .set { ch_mash_fa }

    // Run Mash FA
    MASH_FA (
        ch_mash_fa.fasta,
        ch_mash_fa.fastas
    )

    // BWA channel
    // Contains the reads, contigs, and filtered contigs
    // for each sample combined with each reference
    // selected by Mash FA or the manually-set reference.
    TRIMMOMATIC.out.reads
        .join(SPADES.out.contigs)
        .join(SPADES.out.filtered_contigs)
        .combine(
            MASH_FA.out.fastas
                .map {
                    meta, fastas ->
                    [ meta, meta.set_ref != '' ? [ meta.set_ref ] : fastas ]
                }
                .transpose(),
            by: 0)
        .map {
            meta, reads, contigs, filtered_contigs, fasta ->
            [ meta + [ref: fasta], reads, contigs, filtered_contigs ]
        }
        .cross(
            CHECK_INPUT.out.fasta
                .join(CHECK_INPUT.out.bwa)
                .join(CHECK_INPUT.out.fai)
        ) { it[0].ref }
        .map { it[0] + it[1][1..-1] }
        .multiMap {
            meta, reads, contigs, filtered_contigs, fasta, bwa, fai ->
            reads:            [ meta, reads                                                                  ]
            contigs:          [ meta, contigs                                                                ]
            filtered_contigs: [ meta, filtered_contigs                                                       ]
            fasta:            [ meta, fasta                                                                  ]
            bwa:              [ meta, index                                                                  ]
            fai:              [ meta, fai                                                                    ]
            mapped_threshold: [ meta, meta.set_ref != '' ? 0 : meta.make_ref ? 100 : params.mapped_threshold ]
        }
        .set { ch_bwa }

    // Run BWA
    BWA (
        ch_bwa.reads,
        ch_bwa.fasta,
        ch_bwa.bwa,
        ch_bwa.fai,
        ch_bwa.mapped_threshold
    )

    // BWA output channel
    // Contains the BWA output files and reference files for
    // each sample and its best reference (maximum percent mapped).
    BWA.out.percent_mapped
        .map {
            meta, percent_mapped ->
            [ [id: meta.id], [ meta, percent_mapped ] ]
        }
        .groupTuple()
        .map {
            meta, output ->
            output.max { it[1] }
        }
        .join(BWA.out.depth)
        .join(BWA.out.bam)
        .join(BWA.out.mpileup)
        .join(ch_bwa.contigs)
        .join(ch_bwa.filtered_contigs)
        .join(ch_bwa.fasta)
        .join(ch_bwa.fai)
        .join(ch_bwa.mapped_threshold)
        .multiMap {
            meta, percent_mapped, depth, bam, mpileup, contigs, filtered_contigs, fasta, fai, mapped_threshold ->
            percent_mapped:     [ meta, percent_mapped                                                                         ]
            depth:              [ meta, depth                                                                                  ]
            bam:                [ meta, bam                                                                                    ]
            mpileup:            [ meta, mpileup                                                                                ]
            contigs:            [ meta, contigs                                                                                ]
            filtered_contigs:   [ meta, filtered_contigs                                                                       ]
            fasta:              [ meta, fasta                                                                                  ]
            fai:                [ meta, fai                                                                                    ]
            mapped_threshold:   [ meta, mapped_threshold                                                                       ]
            max_no_ns:          [ meta, meta.set_ref != '' ? 999999999 : meta.make_ref ? 999999999 : params.max_no_ns          ]
            max_no_gaps:        [ meta, meta.set_ref != '' ? 999999999 : meta.make_ref ? 999999999 : params.max_no_gaps        ]
            snp_threshold:      [ meta, meta.set_ref != '' ? 999999999 : meta.make_ref ? 0         : params.snp_threshold      ]
            min_percent_mapped: [ meta, meta.set_ref != '' ? 0         : meta.make_ref ? 0         : params.min_percent_mapped ]
        }
        .set { ch_bwa_out }

    // Run Quast
    QUAST (
        ch_bwa_out.contigs,
        ch_bwa_out.fasta,
        ch_bwa_out.depth,
        ch_bwa_out.percent_mapped,
        ch_bwa_out.max_no_ns,
        ch_bwa_out.max_no_gaps,
        ch_bwa_out.min_percent_mapped,
        ch_bwa_out.mapped_threshold
    )

    // Run Qualimap
    QUALIMAP (
        ch_bwa_out.bam
    )

    // Run Freebayes
    FREEBAYES (
        ch_bwa_out.bam,
        ch_bwa_out.fasta,
        ch_bwa_out.fai,
        ch_bwa_out.snp_threshold,
        QUAST.out.depth_mean,
        QUAST.out.depth_sd
    )

    // Freebayes output channel
    // Branches the samples based on whether they are
    // determined to be close to or distant from the reference
    FREEBAYES.out.mutations
        .join(ch_bwa_out.percent_mapped)
        .join(ch_bwa_out.snp_threshold)
        .join(ch_bwa_out.mapped_threshold)
        .branch {
            meta, mutations, percent_mapped, snp_threshold, mapped_threshold ->
            close: (mutations[0] < snp_threshold) && (percent_mapped > mapped_threshold)
            distant: true
        }
        .set { ch_freebayes_out }

    // Freebayes close channel
    ch_freebayes_out.close
        .join(ch_bwa_out.filtered_contigs)
        .join(ch_bwa_out.fasta)
        .join(ch_bwa_out.mpileup)
        .join(FREEBAYES.out.vcf)
        .multiMap {
            meta, mutations, percent_mapped, snp_threshold, mapped_threshold, filtered_contigs, fasta, mpileup, freebayes ->
            filtered_contigs: [ meta, filtered_contigs ]
            fasta:            [ meta, fasta            ]
            mpileup:          [ meta, mpileup          ]
            freebayes:        [ meta, freebayes        ]
        }
        .set { ch_freebayes_close }

    // Create SNP consensus
    CREATE_SNP_CONS (
        ch_freebayes_close.fasta
            .join(ch_freebayes_close.mpileup)
            .join(ch_freebayes_close.freebayes),
        false
    )

    // Compare SNPs channel
    // Contains the SNP consensus for each sample
    // combined with the list of SNP consensuses
    // of all the samples in its cluster
    CHECK_INPUT.out.snp_cons
        .mix(CHECK_INPUT.out.cluster_snp_cons)
        .groupTuple()
        .cross(
            CREATE_SNP_CONS.out.snp_cons
        ) { it[0].ref }
        .groupTuple()
        .map { [ it[1], it[1].collect { it[1] } + it[0][1] ] }
        .transpose(by: 0)
        .map { it[0] + [ it[1] - it[0][1..-1] ] }
        .multiMap {
            meta, snp_cons, cluster_snp_cons ->
            snp_cons:         [ meta, snp_cons         ]
            cluster_snp_cons: [ meta, cluster_snp_cons ]
        }
        .set { ch_compare_snps }

    // Compare SNPs
    COMPARE_SNPS (
        ch_compare_snps.snp_cons
            .join(ch_compare_snps.cluster_snp_cons)
    )

    // Make MST channel
    // Contains the mutations matrix and
    // list of pairwise diffs for each cluster
    COMPARE_SNPS.out.pairwise_diffs
        .map {
            meta, pairwise_diffs ->
            [ [ref: meta.ref], pairwise_diffs ]
        }
        .groupTuple()
        .join(CHECK_INPUT.out.mutations_matrix)
        .multiMap {
            meta, cluster_pairwise_diffs, mutations_matrix ->
            cluster_pairwise_diffs: [ meta, cluster_pairwise_diffs ]
            mutations_matrix:       [ meta, mutations_matrix       ]
        }
        .set { ch_make_mst }

    // Make MST
    MAKE_MST (
        ch_make_mst.cluster_pairwise_diffs,
        ch_make_mst.mutations_matrix
    )

    // Freebayes distant channel
    ch_freebayes_out.distant
        .join(ch_bwa_out.contigs)
        .join(ch_bwa_out.filtered_contigs)
        .multiMap {
            meta, mutations, percent_mapped, snp_threshold, mapped_threshold, contigs, filtered_contigs ->
            contigs:          [ meta, contigs          ]
            filtered_contigs: [ meta, filtered_contigs ]
        }
        .set { ch_freebayes_distant }

    // Check reference quality
    CHECK_REF_QUAL (
        ch_freebayes_distant.contigs,
        params.med_genome_len
    )

    // Check reference quality output channel
    // Branches the samples based on whether they
    // passed the quality check
    CHECK_REF_QUAL.out.csv
        .map {
            meta, csv ->
            [ meta ] + csv.splitCsv().collect()
        }
        .branch {
            meta, passed_qc ->
            passed_qc: passed_qc[0].toBoolean()
            failed_qc: true
        }
        .set { ch_check_ref_qual_out }

    // Make reference channel
    // Contains the assembled fasta for each sample
    // that was determined to be distant
    // and passed the quality check
    ch_check_ref_qual_out.passed_qc
        .join(ch_freebayes_distant.filtered_contigs)
        .map {
            meta, passed_qc, fasta ->
            [ [ref: meta.ref], fasta ]
        }
        .set { ch_make_reference }

    // Make reference
    MAKE_REFERENCE (
        ch_make_reference
    )

    // Make references reference channel
    // Contains information for each reference
    MAKE_REFERENCE.out.fasta
        .mix(CHECK_INPUT.out.fasta)
        .join(
            MAKE_REFERENCE.out.snp_cons
                .mix(CHECK_INPUT.out.snp_cons)
        )
        .join(
            MAKE_REFERENCE.out.bwa
                .mix(CHECK_INPUT.out.bwa)
        )
        .join(
            MAKE_REFERENCE.out.fai
                .mix(CHECK_INPUT.out.fai)
        )
        .join(
            MAKE_REFERENCE.out.mutations_matrix
                .mix(
                    MAKE_MST.out.mutations_matrix
                        .concat(CHECK_INPUT.out.mutations_matrix)
                        .unique { it[0].ref }
                )
        )
        .map {
            meta, fasta, snp_cons, bwa, fai, mutations_matrix ->
            [ meta.ref, meta.ref, fasta, snp_cons, bwa, fai, mutations_matrix ]
        }
        .set { ch_make_references_reference }

    // Make references cluster reference channel
    // Contains information for each cluster reference
    ch_freebayes_close.filtered_contigs
        .mix(CHECK_INPUT.out.cluster_fasta)
        .join(
            CREATE_SNP_CONS.out.snp_cons
                .mix(CHECK_INPUT.out.cluster_snp_cons)
        )
        .map {
            meta, fasta, snp_cons ->
            [ meta.id, meta.ref, fasta, snp_cons, '', '', '' ]
        }
        .set { ch_make_references_cluster_reference }

    // Make references
    MAKE_REFERENCES (
        ch_make_references_reference
            .mix(ch_make_references_cluster_reference)
            .groupTuple(by: [])
            .dump(tag: 'make_references')
    )

    // Collect reports
    ch_reports = ch_reports.concat(TRIMMOMATIC.out.reports)
    ch_reports = ch_reports.concat(FASTQC.out.reports)
    ch_reports = ch_reports.concat(MASH_FQ.out.reports)
    ch_reports = ch_reports.concat(SPADES.out.reports)
    ch_reports = ch_reports.concat(MASH_FA.out.reports)
    ch_reports = ch_reports.concat(BWA.out.reports)
    ch_reports = ch_reports.concat(QUAST.out.reports)
    ch_reports = ch_reports.concat(QUALIMAP.out.reports)
    ch_reports = ch_reports.concat(FREEBAYES.out.reports)

    // Create report
    CREATE_REPORT (
        ch_reports
            .map {
                meta, report ->
                [ meta - [ref: meta.ref], report ]
            }
            .groupTuple()
            .join(CHECK_INPUT.out.reads),
        params.genome
    )

    // Collect versions
    ch_versions = ch_versions.mix(CHECK_INPUT.out.versions)
    ch_versions = ch_versions.mix(TRIMMOMATIC.out.versions)
    ch_versions = ch_versions.mix(FASTQC.out.versions)
    ch_versions = ch_versions.mix(MASH_FQ.out.versions)
    ch_versions = ch_versions.mix(SPADES.out.versions)
    ch_versions = ch_versions.mix(MASH_FA.out.versions)
    ch_versions = ch_versions.mix(BWA.out.versions)
    ch_versions = ch_versions.mix(QUAST.out.versions)
    ch_versions = ch_versions.mix(QUALIMAP.out.versions)
    ch_versions = ch_versions.mix(FREEBAYES.out.versions)

    // Dump versions
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
