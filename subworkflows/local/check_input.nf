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
        .set { reads }

    CHECK_REFERENCES (
        references,
        true
    )

    CHECK_REFERENCES.out.csv
        .splitCsv(header: true, sep: ',')
        .map { create_fasta_channel(it) }
        .set { fasta }

    // Collect versions
    ch_versions = ch_versions.mix(CHECK_SAMPLES.out.versions)
    ch_versions = ch_versions.mix(CHECK_REFERENCES.out.versions)

    emit:
    reads                  // channel: [ val(meta), [ reads ] ]
    fasta                  // channel: [ val(meta), [ fasta ] ]
    versions = ch_versions // channel: [ versions.yml ]
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
def create_fasta_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.ref = row.reference

    if (!file(row.fasta).exists()) {
        exit 1, "ERROR: Please check reference samplesheet -> FASTA file does not exist!\n${row.fasta}"
    }
    return [ meta, file(row.fasta) ]
}
