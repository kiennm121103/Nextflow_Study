// This files contains the main Nextflow pipeline definition

// GENOME PARAMETERS VALUES
params.fasta = getGenomeAttribute('fasta')

// IMPORT FUNCTIONS/ MODULES/ SUBWORKFLOWS/ WORKFLOWS

include {RNASEQ} from "./workflows/rnaseq"
include {}


// FUNCTIONS

// Get attribute from genome config file e.g. fasta

def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}