//
// This file holds several functions specific to the workflow/generate_references.nf in the nf-core/legiocluster pipeline
//

class WorkflowGenerateReferences {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {
        if (!params.references) {
            log.error "References path not specified!"
            System.exit(1)
        }
    }

    //
    // Get genomes FASTA list
    //
    public static List genomesFastaList(params) {
        List genome_fastas = []
        for (genome in params.genomes.keySet()) {
            if (params.genomes.get(genome).containsKey('fasta')) {
                genome_fastas.add(params.genomes.get(genome).fasta)
            }
        }
        return genome_fastas
    }
}
