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
            reference: it.sample == it.reference
            cluster_reference: true
        }
        .set { ch_reference_branch }

    ch_reference_branch.reference
        .map { create_reference_channel(it, false) }
        .multiMap {
            meta, fasta, snp_cons, bwa, fai, mutations_matrix ->
            fasta:            [ meta, fasta            ]
            snp_cons:         [ meta, snp_cons         ]
            bwa:              [ meta, bwa              ]
            fai:              [ meta, fai              ]
            mutations_matrix: [ meta, mutations_matrix ]
        }
        .set { ch_reference }

    ch_reference_branch.cluster_reference
        .map { create_reference_channel(it, true) }
        .multiMap {
            meta, fasta, snp_cons ->
            fasta:    [ meta, fasta    ]
            snp_cons: [ meta, snp_cons ]
        }
        .set { ch_cluster_reference }

    ch_reference.fasta
        .map {
            meta, fasta ->
            [ [id: meta.ref], fasta ]
        }
        .mix(ch_cluster_reference.fasta)
        .cross(ch_reads) { it[0].id }
        .map { error "[ERROR] An isolate with the name '${it[0][0].id}' already exists as a reference." }

    ch_reads
        .count()
        .filter { it > 1 }
        .view { "[NOTICE] Multiple samples will cluster independently of each other. If this behavior is not desired, run each sample sequentially instead." }

    // Collect versions
    ch_versions = ch_versions.mix(CHECK_SAMPLES.out.versions)
    ch_versions = ch_versions.mix(CHECK_REFERENCES.out.versions)

    emit:
    reads            = ch_reads
    fasta            = ch_reference.fasta
    snp_cons         = ch_reference.snp_cons
    bwa              = ch_reference.bwa
    fai              = ch_reference.fai
    mutations_matrix = ch_reference.mutations_matrix
    cluster_fasta    = ch_cluster_reference.fasta
    cluster_snp_cons = ch_cluster_reference.snp_cons
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
def create_reference_channel(LinkedHashMap row, boolean cluster_reference) {
    // create meta map
    def meta = [:]
    meta.ref = row.reference

    if (!file(row.fasta).exists()) {
        exit 1, "ERROR: Please check reference samplesheet -> FASTA file does not exist!\n${row.fasta}"
    }
    if (!file(row.snp_cons).exists()) {
        exit 1, "ERROR: Please check reference samplesheet -> SNP consensus file does not exist!\n${row.fasta}"
    }

    if (!cluster_reference) {
        if (!file(row.bwa, type: 'dir').exists()) {
            exit 1, "ERROR: Please check reference samplesheet -> BWA directory does not exist!\n${row.fasta}"
        }
        if (!file(row.fai).exists()) {
            exit 1, "ERROR: Please check reference samplesheet -> FAI file does not exist!\n${row.fasta}"
        }
        if (!file(row.mutations_matrix).exists()) {
            exit 1, "ERROR: Please check reference samplesheet -> Mutations matrix file does not exist!\n${row.fasta}"
        }
        return [ meta, file(row.fasta), file(row.snp_cons), file(row.bwa, type: 'dir'), file(row.fai), file(row.mutations_matrix) ]
    }

    meta.id = row.sample
    return [ meta, file(row.fasta), file(row.snp_cons) ]
}
