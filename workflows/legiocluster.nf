/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowLegiocluster.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

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
include { MASH_SKETCH as MASH_SKETCH_SPECIES } from '../modules/local/mash_sketch'
include { MASH_SKETCH as MASH_SKETCH_STRAINS } from '../modules/local/mash_sketch'
include { CREATE_REPORT                      } from '../modules/local/create_report'
include { CUSTOM_DUMPSOFTWAREVERSIONS        } from '../modules/local/custom_dumpsoftwareversions'
include { MULTIQC                            } from '../modules/local/multiqc'

// Subworkflows
include { CHECK_INPUT     } from '../subworkflows/local/check_input'
include { RUN_TRIMMOMATIC } from '../subworkflows/local/run_trimmomatic'
include { RUN_FASTQC      } from '../subworkflows/local/run_fastqc'
include { RUN_MASH_FQ     } from '../subworkflows/local/run_mash_fq'
include { RUN_SPADES      } from '../subworkflows/local/run_spades'
include { RUN_MASH_FA     } from '../subworkflows/local/run_mash_fa'
include { RUN_BWA_FA      } from '../subworkflows/local/run_bwa_fa'
include { RUN_BWA         } from '../subworkflows/local/run_bwa'

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
        ch_input
    )

    // Run trimmomatic
    RUN_TRIMMOMATIC (
        CHECK_INPUT.out.reads
    )

    // Run fastqc
    RUN_FASTQC (
        CHECK_INPUT.out.reads,
        RUN_TRIMMOMATIC.out.reads
    )

    // Mash sketch (species)
    MASH_SKETCH_SPECIES (
        Channel.fromPath(params.species_refs).collect().map { [ [:], it ] }
    )

    // Run mash fq
    RUN_MASH_FQ (
        RUN_TRIMMOMATIC.out.reads,
        MASH_SKETCH_SPECIES.out.mash.map { it[1] }
    )

    // Run spades
    RUN_SPADES (
        RUN_TRIMMOMATIC.out.reads,
        RUN_TRIMMOMATIC.out.max_read_len
    )

    // Create strain fasta channel [ meta(ref), fasta ]
    Channel.fromPath(params.strain_refs)
        .map {
            fasta ->
            [ [ref: fasta.name], fasta ]
        }
        .set { ch_strain_fasta }

    // Mash sketch (strains)
    MASH_SKETCH_STRAINS (
        ch_strain_fasta
            .collect { it[1] }
            .map { [ [:], it ] }
    )

    // Run mash fa
    RUN_MASH_FA (
        RUN_TRIMMOMATIC.out.reads,
        RUN_SPADES.out.fasta,
        MASH_SKETCH_STRAINS.out.mash.map { it[1] }
    )

    // Create bwa reads, fasta channel [ meta(id, ref), reads, fasta ]
    RUN_TRIMMOMATIC.out.reads
        .combine(
            RUN_MASH_FA.out.fastas
                .map {
                    meta, fastas ->
                    [ meta, meta.set_ref != 'NO_FILE' ? [ meta.set_ref ] : fastas ]
                }
                .transpose(),
            by: 0)
        .map {
            meta, reads, fasta ->
            [ meta + [ref: fasta], reads ]
        }
        .cross(ch_strain_fasta) { it[0].ref }
        .map { it[0] + [ it[1][1] ] }
        .set { ch_bwa_reads_fasta }

    // Run bwa fa
    RUN_BWA_FA (
        ch_bwa_reads_fasta
            .map {
                meta, reads, fasta ->
                [ [ref: meta.ref], fasta ]
            }
            .unique()
    )

    // Create bwa reads, fasta, index, fai channels
    ch_bwa_reads_fasta
        .cross(
            RUN_BWA_FA.out.index
                .join(RUN_BWA_FA.out.fai)
        ) { it[0].ref }
        .map { it[0] + [ it[1][1], it[1][2] ] }
        .multiMap {
            meta, reads, fasta, index, fai ->
            reads: [ meta, reads ]
            fasta: [ meta, fasta ]
            index: [ meta, index ]
            fai: [ meta, fai ]
        }
        .set { ch_bwa }

    // Run bwa
    RUN_BWA (
        ch_bwa.reads,
        ch_bwa.fasta,
        ch_bwa.index,
        ch_bwa.fai
    )

    // Collect reports
    ch_reports = ch_reports.concat(RUN_TRIMMOMATIC.out.reports)
    ch_reports = ch_reports.concat(RUN_FASTQC.out.reports)
    ch_reports = ch_reports.concat(RUN_MASH_FQ.out.reports)
    ch_reports = ch_reports.concat(RUN_SPADES.out.reports)
    ch_reports = ch_reports.concat(RUN_MASH_FA.out.reports)
    ch_reports = ch_reports.concat(RUN_BWA.out.reports)

    CREATE_REPORT (
        ch_reports.groupTuple().join(CHECK_INPUT.out.reads)
    )

    // Collect versions
    ch_versions = ch_versions.mix(CHECK_INPUT.out.versions)
    ch_versions = ch_versions.mix(RUN_TRIMMOMATIC.out.versions)
    ch_versions = ch_versions.mix(RUN_FASTQC.out.versions)
    ch_versions = ch_versions.mix(RUN_MASH_FQ.out.versions)
    ch_versions = ch_versions.mix(RUN_SPADES.out.versions)
    ch_versions = ch_versions.mix(RUN_MASH_FA.out.versions)
    ch_versions = ch_versions.mix(RUN_BWA_FA.out.versions)
    ch_versions = ch_versions.mix(RUN_BWA.out.versions)

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
    ch_multiqc_files = ch_multiqc_files.mix(RUN_FASTQC.out.zip.collect{ it[1] }.ifEmpty([]))

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
