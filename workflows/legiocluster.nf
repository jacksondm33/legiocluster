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
include { CHECK_INPUT } from '../subworkflows/local/check_input'
include { TRIMMOMATIC } from '../subworkflows/local/trimmomatic'
include { FASTQC      } from '../subworkflows/local/fastqc'
include { MASH_FQ     } from '../subworkflows/local/mash_fq'
include { SPADES      } from '../subworkflows/local/spades'
include { MASH_FA     } from '../subworkflows/local/mash_fa'
include { BWA_FA      } from '../subworkflows/local/bwa_fa'
include { BWA         } from '../subworkflows/local/bwa'
include { QUAST       } from '../subworkflows/local/quast'

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
    TRIMMOMATIC (
        CHECK_INPUT.out.reads
    )

    // Run fastqc
    FASTQC (
        CHECK_INPUT.out.reads,
        TRIMMOMATIC.out.reads
    )

    // Mash sketch (species)
    MASH_SKETCH_SPECIES (
        Channel.fromPath(params.species_refs).collect().map { [ [:], it ] }
    )

    // Run mash fq
    MASH_FQ (
        TRIMMOMATIC.out.reads,
        MASH_SKETCH_SPECIES.out.mash.map { it[1] }
    )

    // Run spades
    SPADES (
        TRIMMOMATIC.out.reads,
        TRIMMOMATIC.out.max_read_len
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
    MASH_FA (
        TRIMMOMATIC.out.reads,
        SPADES.out.contigs,
        MASH_SKETCH_STRAINS.out.mash.map { it[1] }
    )

    // Create bwa reads, fasta channel [ meta(id, ref), reads, fasta ]
    TRIMMOMATIC.out.reads
        .combine(
            MASH_FA.out.fastas
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
    BWA_FA (
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
        }
        .set { ch_bwa }

    // Run bwa
    BWA (
        ch_bwa.reads,
        ch_bwa.fasta,
        ch_bwa.index,
        ch_bwa.fai
    )

    BWA.out.percent_mapped
        .join(BWA.out.depth)
        .cross(ch_strain_fasta) { it[0].ref }
        .map { it[0] + [ it[1][1] ] }
        .map {
            meta, percent_mapped, depth, fasta ->
            [ meta - [ref: meta.ref], [percent_mapped, depth, fasta] ]
        }
        .groupTuple()
        .map {
            meta, percent_mapped_depth_fasta ->
            [ meta ] + percent_mapped_depth_fasta.max { it[0] }
        }
        .multiMap {
            meta, percent_mapped, depth, fasta ->
            percent_mapped: [ meta, percent_mapped ]
            depth: [ meta, depth ]
            fasta: [ meta, fasta ]
        }
        .set { ch_quast }

    // Run quast
    QUAST (
        SPADES.out.contigs,
        ch_quast.fasta,
        ch_quast.depth,
        ch_quast.percent_mapped
    )

    // Collect reports
    ch_reports = ch_reports.concat(TRIMMOMATIC.out.reports)
    ch_reports = ch_reports.concat(FASTQC.out.reports)
    ch_reports = ch_reports.concat(MASH_FQ.out.reports)
    ch_reports = ch_reports.concat(SPADES.out.reports)
    ch_reports = ch_reports.concat(MASH_FA.out.reports)
    ch_reports = ch_reports.concat(BWA.out.reports)
    ch_reports = ch_reports.concat(QUAST.out.reports)

    CREATE_REPORT (
        ch_reports
            .map {
                meta, report ->
                [ meta - [ref: meta.ref], report ]
            }
            .groupTuple()
            .join(CHECK_INPUT.out.reads)
    )

    // Collect versions
    ch_versions = ch_versions.mix(CHECK_INPUT.out.versions)
    ch_versions = ch_versions.mix(TRIMMOMATIC.out.versions)
    ch_versions = ch_versions.mix(FASTQC.out.versions)
    ch_versions = ch_versions.mix(MASH_FQ.out.versions)
    ch_versions = ch_versions.mix(SPADES.out.versions)
    ch_versions = ch_versions.mix(MASH_FA.out.versions)
    ch_versions = ch_versions.mix(BWA_FA.out.versions)
    ch_versions = ch_versions.mix(BWA.out.versions)
    ch_versions = ch_versions.mix(QUAST.out.versions)

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
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{ it[1] }.ifEmpty([]))

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
