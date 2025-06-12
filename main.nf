#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    SEQAUDIT - SEQUENCING DATA SIMULATION PIPELINE
========================================================================================
    A comprehensive pipeline to simulate ONT, PacBio, and Illumina sequencing data
    from reference genomes

    Usage:
    nextflow run main.nf --input samplesheet.csv --outdir results

    Input formats supported:
    - RefSeq ID (e.g., GCF_000005825.2)
    - GenBank ID (e.g., GCA_000001405.29)
    - Local FASTA file path
----------------------------------------------------------------------------------------
*/

// Include main subworkflow
include { SIMULATE_READS } from './subworkflows/local/simulate_reads'

// Print help message if requested
if (params.help) {
    log.info """
    Usage:
    nextflow run main.nf --input samplesheet.csv --outdir results

    Required arguments:
      --input                   Path to comma-separated file containing information about the samples
      --outdir                  The output directory where the results will be saved

    Optional arguments:
      --ont_simulator           ONT simulator to use ('pbsim3' or 'nanosim') [default: pbsim3]
      --pacbio_simulator        PacBio simulator to use ('pbsim3') [default: pbsim3]
      --help                    Show this help message and exit
    """
    exit 0
}

// Print parameter summary
log.info """
=======================================================
SEQAUDIT PIPELINE PARAMETERS
=======================================================
Input/Output:
  input                     : ${params.input}
  outdir                    : ${params.outdir}

ONT Simulation:
  ont_simulator             : ${params.ont_simulator}
  ont_model_url             : ${params.ont_model_url}

PacBio Simulation:
  pacbio_simulator          : ${params.pacbio_simulator}
  pacbio_model_url          : ${params.pacbio_model_url}

Illumina Simulation:
  illumina_read_length      : ${params.illumina_read_length}
  illumina_fragment_mean    : ${params.illumina_fragment_mean}
=======================================================
"""

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {

    // Create input channel from samplesheet
    ch_input = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [:]
            meta.id = row.sample_id
            meta.genome_source = row.genome_source // 'refseq', 'genbank', or 'local'
            meta.genome_id = row.genome_id // RefSeq/GenBank ID or file path
            meta.ont_reads = row.ont_reads as Integer ?: 0
            meta.pacbio_reads = row.pacbio_reads as Integer ?: 0
            meta.illumina_reads = row.illumina_reads as Integer ?: 0

            return [meta, row.genome_id]
        }

    // Run simulation subworkflow
    SIMULATE_READS(ch_input)

    // Emit outputs
    emit:
    ont_reads = SIMULATE_READS.out.ont_reads
    pacbio_reads = SIMULATE_READS.out.pacbio_reads
    illumina_reads = SIMULATE_READS.out.illumina_reads
    qc_summary_table = SIMULATE_READS.out.qc_summary_table
}
