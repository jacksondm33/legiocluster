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
include { MASH_SKETCH as MASH_SKETCH_SPECIES } from '../modules/local/mash_sketch'
include { MAKE_SNP_CONS                      } from '../modules/local/make_snp_cons'
include { COMPARE_SNPS                       } from '../modules/local/compare_snps'
include { CHECK_REF_QUAL                     } from '../modules/local/check_ref_qual'
include { CONVERT_REPORTS                    } from '../modules/local/convert_reports'
include { MAKE_SOFTWARE_VERSIONS             } from '../modules/local/make_software_versions'
include { MAKE_REPORT                        } from '../modules/local/make_report'
include { MULTIQC                            } from '../modules/local/multiqc'
include { MULTIQC as MULTIQC_SUMMARY         } from '../modules/local/multiqc'

// Subworkflows
include { CHECK_INPUT    } from '../subworkflows/local/check_input'
include { TRIMMOMATIC    } from '../subworkflows/local/trimmomatic'
include { FASTQC         } from '../subworkflows/local/fastqc'
include { MASH_FQ        } from '../subworkflows/local/mash_fq'
include { SPADES         } from '../subworkflows/local/spades'
include { MASH_FA        } from '../subworkflows/local/mash_fa'
include { BWA            } from '../subworkflows/local/bwa'
include { KRAKEN         } from '../subworkflows/local/kraken'
include { QUAST          } from '../subworkflows/local/quast'
include { QUALIMAP       } from '../subworkflows/local/qualimap'
include { FREEBAYES      } from '../subworkflows/local/freebayes'
include { MAKE_MST       } from '../subworkflows/local/make_mst'
include { PARSNP         } from '../subworkflows/local/parsnp'
include { MAKE_REFERENCE } from '../subworkflows/local/make_reference'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GENERATE_REFERENCES {

    // Mash sketch (species)
    MASH_SKETCH_SPECIES (
        ch_genomes.map { [ [:], it ] }
    )

    ch_references
        .map { [ [ref: it.baseName], it ] }
        .set { ch_fasta }

    MAKE_REFERENCE (
        ch_fasta
    )

    MAKE_REFERENCE.out.fasta
        .join(MAKE_REFERENCE.out.snp_cons)
        .join(MAKE_REFERENCE.out.bwa)
        .join(MAKE_REFERENCE.out.fai)
        .join(MAKE_REFERENCE.out.mutations_matrix)
        .map {
            meta, fasta, snp_cons, bwa, fai, mutations_matrix ->
            [ meta.ref, meta.ref, fasta, snp_cons, bwa, fai, mutations_matrix ].join(',')
        }
        .set { ch_make_references }

    references_header = [ 'sample', 'reference', 'fasta', 'snp_cons', 'bwa', 'fai', 'mutations_matrix' ].join(',')
    ch_make_references.collectFile(name: "references_${params.genome}.csv", newLine: true, seed: references_header, sort: true, storeDir: params.outdir)

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
