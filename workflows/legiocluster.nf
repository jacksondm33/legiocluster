/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowLegiocluster.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.references, params.illuminaclip, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Create input and references channels
ch_input = file(params.input)
ch_references = file(params.references)

// Create genomes channel
ch_genomes = Channel.value(WorkflowLegiocluster.genomesFastaList(params))
    .map { it.collect { file(it, checkIfExists: true) } }

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

include { LEGIOCLUSTER_MAIN } from './legiocluster_main'

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

    // Mash sketch
    MASH_SKETCH (
        ch_genomes.map { [ [:], it ] },
        false,
        true
    )

    // Run mash fq
    MASH_FQ (
        TRIMMOMATIC.out.reads,
        MASH_SKETCH.out.mash.map { it[1] }
    )

    // Run spades
    SPADES (
        TRIMMOMATIC.out.reads,
        TRIMMOMATIC.out.max_read_len
    )

    // Hack to make nextflow recursion work
    TRIMMOMATIC.out.reads
        .combine(
            CHECK_INPUT.out.fasta
                .collect { it[1] }
                .map { [ it ] }
        )
        .map {
            meta, reads, fastas ->
            [ meta, fastas ]
        }
        .set { ch_fastas }

    LEGIOCLUSTER_MAIN
        .scan(TRIMMOMATIC.out.reads, SPADES.out.contigs, SPADES.out.filtered_contigs, CHECK_INPUT.out.fasta, Channel.empty(), Channel.empty())

    // Create SNP consensus (fasta)
    CREATE_SNP_CONS_FA (
        ch_freebayes_close.fasta
            .map {
                meta, fasta ->
                [ [ref: meta.ref], fasta ]
            }
            .unique()
    )

    // Create SNP consensus input channel
    ch_freebayes_close.mpileup
        .join(ch_freebayes_close.freebayes)
        .cross(
            CREATE_SNP_CONS_FA.out.bases
        ) { it[0].ref }
        .map { it[0] + it[1][1..-1] }
        .set { ch_create_snp_cons_in }

    // Create SNP consensus
    CREATE_SNP_CONS (
        ch_create_snp_cons_in,
        false
    )

    // Compare SNPs input channel
    CREATE_SNP_CONS_FA.out.snp_cons
        .cross(
            CREATE_SNP_CONS.out.snp_cons
        ) { it[0].ref }
        .groupTuple()
        .map { [ it[1], it[1].collect { it[1] } + [ it[0][1] ] ] }
        .transpose(by: 0)
        .map { it[0] + [ it[1] - [ it[0][1] ] ] }
        .set { ch_compare_snps_in }

    COMPARE_SNPS (
        ch_compare_snps_in
    )

    COMPARE_SNPS.out.pairwise_diffs
        .map {
            meta, pairwise_diffs ->
            [ [ref: meta.ref], pairwise_diffs ]
        }
        .groupTuple()
        .set { ch_make_mst_in }

    MAKE_MST (
        ch_make_mst_in
    )

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
