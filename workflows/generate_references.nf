/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Validate input parameters
WorkflowGenerateReferences.initialise(params, log)

// Create input channels
ch_references = Channel.fromPath(params.references + '/*.fa', checkIfExists: true)

// Create genomes channel
ch_genomes = Channel.value(WorkflowGenerateReferences.genomesFastaList(params))
    .map { it.collect { file(it, checkIfExists: true) } }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Modules
include { MASH_SKETCH    } from '../modules/local/mash_sketch'

// Subworkflows
include { MAKE_REFERENCE } from '../subworkflows/local/make_reference'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GENERATE_REFERENCES {

    // Mash sketch (genomes)
    MASH_SKETCH (
        ch_genomes.map { [ [:], it ] },
        false,
        true
    )

    ch_references
        .map { [ [ref: it.baseName], it ] }
        .set { ch_fasta }

    MAKE_REFERENCE (
        ch_fasta
    )

    MAKE_REFERENCE.out.fasta
        .join(MAKE_REFERENCE.out.bwa)
        .join(MAKE_REFERENCE.out.fai)
        .join(MAKE_REFERENCE.out.snp_cons)
        .map {
            meta, fasta, bwa, fai, snp_cons ->
            meta.ref + ',' + fasta + ',' + bwa + ',' + fai + ',' + snp_cons
        }
        .collectFile(name: "references_${params.genome}.csv", newLine: true, sort: true, storeDir: params.outdir)

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
