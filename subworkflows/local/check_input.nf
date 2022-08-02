//
// Check input samplesheet and get read channels
//

include { CHECK_SAMPLESHEET } from '../../modules/local/check_samplesheet'

workflow CHECK_INPUT {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    CHECK_SAMPLESHEET ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it) }
        .set { reads }

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
    versions = CHECK_SAMPLESHEET.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    meta.set_ref    = row.set_ref
    meta.make_ref   = row.make_ref.toBoolean()

    // add path(s) of the fastq file(s) to the meta map
    if (row.fastq_1 == '') {
        return [ meta, [] ]
    }
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (!file(row.fastq_2).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
    }
    return [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
}
