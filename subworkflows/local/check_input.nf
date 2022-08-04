//
// Check input samplesheet and get read channels
//

include { CHECK_INPUT as CHECK_SAMPLES    } from '../../modules/local/check_input'
include { CHECK_INPUT as CHECK_REFERENCES } from '../../modules/local/check_input'

workflow CHECK_INPUT {
    take:
    samples
    references

    main:
    ch_versions = Channel.empty()

    CHECK_SAMPLES (
        samples,
        false
    )

    CHECK_SAMPLES.out.csv
        .splitCsv(header: true, sep: ',')
        .map { create_reads_channel(it) }
        .set { ch_reads }

    CHECK_REFERENCES (
        references,
        true
    )

    CHECK_REFERENCES.out.csv
        .splitCsv(header: true, sep: ',')
        .branch {
            reference: it.reference == it.cluster
            cluster: true
        }
        .set { ch_reference }

    ch_reference.reference
        .map { create_reference_channel(it) }
        .multiMap {
            meta, fasta, bwa, fai, snp_cons, mutations_matrix ->
            fasta:            [ meta, fasta            ]
            bwa:              [ meta, bwa              ]
            fai:              [ meta, fai              ]
            snp_cons:         [ meta, snp_cons         ]
            mutations_matrix: [ meta, mutations_matrix ]
        }
        .set { ch_output }

    ch_reference.cluster
        .map { create_cluster_channel(it) }
        .multiMap {
            meta, fasta, snp_cons ->
            fasta:    [ meta, fasta    ]
            snp_cons: [ meta, snp_cons ]
        }
        .set { ch_cluster_output }

    // Collect versions
    ch_versions = ch_versions.mix(CHECK_SAMPLES.out.versions)
    ch_versions = ch_versions.mix(CHECK_REFERENCES.out.versions)

    emit:
    reads            = ch_reads
    fasta            = ch_output.fasta
    bwa              = ch_output.bwa
    fai              = ch_output.fai
    snp_cons         = ch_output.snp_cons
    mutations_matrix = ch_output.mutations_matrix
    cluster_fasta    = ch_cluster_output.fasta
    cluster_snp_cons = ch_cluster_output.snp_cons
    versions         = ch_versions                // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_reads_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    meta.set_ref    = row.set_ref
    meta.make_ref   = row.make_ref.toBoolean()

    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FASTQ file does not exist!\n${row.fastq_1}"
    }
    if (!file(row.fastq_2).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 2 FASTQ file does not exist!\n${row.fastq_2}"
    }
    return [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
}

// Function to get list of [ meta, fasta ]
def create_reference_channel(LinkedHashMap row, boolean cluster) {
    // create meta map
    def meta = [:]
    meta.id  = row.reference
    meta.ref = row.cluster

    if (!file(row.fasta).exists()) {
        exit 1, "ERROR: Please check reference samplesheet -> FASTA file does not exist!\n${row.fasta}"
    }
    if (!file(row.snp_cons).exists()) {
        exit 1, "ERROR: Please check reference samplesheet -> SNP consensus file does not exist!\n${row.fasta}"
    }
    if (!cluster) {
        if (!file(row.bwa, type: 'dir').exists()) {
            exit 1, "ERROR: Please check reference samplesheet -> BWA directory does not exist!\n${row.fasta}"
        }
        if (!file(row.fai).exists()) {
            exit 1, "ERROR: Please check reference samplesheet -> FAI file does not exist!\n${row.fasta}"
        }
        if (!file(row.mutations_matrix).exists()) {
            exit 1, "ERROR: Please check reference samplesheet -> Mutations matrix file does not exist!\n${row.fasta}"
        }
        return [ meta, file(row.fasta), file(row.bwa, type: 'dir'), file(row.fai), file(row.snp_cons), file(row.mutations_matrix) ]
    }
    return [ meta, file(row.fasta), file(row.snp_cons) ]
}
