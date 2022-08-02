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
include { CREATE_SNP_CONS_FA                 } from '../modules/local/create_snp_cons_fa'
include { CREATE_SNP_CONS                    } from '../modules/local/create_snp_cons'
include { COMPARE_SNPS                       } from '../modules/local/compare_snps'
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

    // Create species fastas channel [ meta(), [ fastas ] ]
    Channel.fromPath(params.species_refs)
        .collect()
        .map { [ [:], it ] }
        .set { ch_species_fastas }

    // Mash sketch (species)
    MASH_SKETCH_SPECIES (
        ch_species_fastas,
        false,
        true
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

    SPADES.out.filtered_contigs
        .branch {
            meta, filtered_contigs ->
            make_ref: meta.make_ref
            set_ref: meta.set_ref != ''
            get_ref: true
        }
        .set { ch_spades_out }

    // Channel.fromPath(params.strain_refs)
    //     .map { [ [ref: it.baseName], it ] }
    //     .set { ch_strain_refs }

    ch_spades_out.make_ref
        .map {
            meta, reads, contigs, filtered_contigs ->
            filtered_contigs
        }
        .mix(Channel.fromPath(params.strain_refs))
        .collect()
        .set { ch_refs }

    ch_spades_out.get_ref
        .multiMap {
            meta, reads, contigs, filtered_contigs ->
            reads: [ meta, reads ]
            contigs: [ meta, contigs ]
            filtered_contigs: [ meta, filtered_contigs ]
        }
        .set { ch_main }

    LEGIOCLUSTER_MAIN
        .scan(ch_main.reads, ch_main.contigs, ch_main.filtered_contigs, ch_refs, Channel.empty(), Channel.empty())

    LEGIOCLUSTER_MAIN.out.fastas.dump(tag: 'final')

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
