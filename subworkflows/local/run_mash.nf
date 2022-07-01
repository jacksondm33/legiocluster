include { MASH_SKETCH_SPECIES } from '../../modules/local/mash_sketch_species'
include { CONCATENATE_READS } from '../../modules/local/concatenate_reads'
include { MASH_SKETCH } from '../../modules/local/mash_sketch'
include { MASH_DIST } from '../../modules/nf-core/modules/mash/dist/main'
include { PARSE_MASH_OUTPUT } from '../../modules/local/parse_mash_output'

workflow RUN_MASH {
    take:
    reads
    report

    main:
    ch_versions = Channel.empty()

    MASH_SKETCH_SPECIES (
        Channel.fromPath(params.species_ref).collect()
    )

    CONCATENATE_READS (
        reads
    )

    MASH_SKETCH (
        CONCATENATE_READS.out.comb_reads
    )

    MASH_DIST (
        MASH_SKETCH.out.mash,
        MASH_SKETCH_SPECIES.out.mash
    )

    PARSE_MASH_OUTPUT (
        MASH_DIST.out.dist.join(report),
        Channel.of(params.species.collect { sp_abbr, info -> sp_abbr + ',' + info[0] }).flatten().collectFile(name: "species.csv", newLine: true, sort: true),
        Channel.of('RvSp')
    )

    // Collect versions
    ch_versions = ch_versions.mix(MASH_SKETCH_SPECIES.out.versions)
    ch_versions = ch_versions.mix(MASH_SKETCH.out.versions)
    ch_versions = ch_versions.mix(MASH_DIST.out.versions)
    ch_versions = ch_versions.mix(PARSE_MASH_OUTPUT.out.versions)

    emit:
    report = report
    versions = ch_versions // channel: [ versions.yml ]
}
