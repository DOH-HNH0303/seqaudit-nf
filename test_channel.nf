#!/usr/bin/env nextflow

// Simple test to verify channel branching works correctly
nextflow.enable.dsl = 2

workflow {
    // Create test data similar to your pipeline
    test_data = Channel.of(
        [[id: 'test_ecoli', ont_reads: 1, pacbio_reads: 1, illumina_reads: 1], 'genome1.fasta'],
        [[id: 'test_small', ont_reads: 0, pacbio_reads: 1, illumina_reads: 1], 'genome2.fasta']
    )

    // Test the branching approach
    test_data
        .branch { meta, genome ->
            ont: meta.ont_reads > 0
                println "ONT: ${meta.id} - reads: ${meta.ont_reads}"
                return true
            pacbio: meta.pacbio_reads > 0
                println "PacBio: ${meta.id} - reads: ${meta.pacbio_reads}"
                return true
            illumina: meta.illumina_reads > 0
                println "Illumina: ${meta.id} - reads: ${meta.illumina_reads}"
                return true
        }
        .set { branched }

    // Count items in each branch
    branched.ont.count().view { "ONT samples: $it" }
    branched.pacbio.count().view { "PacBio samples: $it" }
    branched.illumina.count().view { "Illumina samples: $it" }
}
