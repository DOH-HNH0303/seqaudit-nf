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

// Include nf-core utility functions
include { UTILS_NEXTFLOW_PIPELINE } from './subworkflows/nf-core/utils_nextflow_pipeline'
include { UTILS_NFCORE_PIPELINE } from './subworkflows/nf-core/utils_nfcore_pipeline'

// Include main subworkflow
include { SIMULATE_READS } from './subworkflows/local/simulate_reads'

/*
========================================================================================
    NAMED WORKFLOWS FOR PIPELINE
========================================================================================
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow SEQAUDIT {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()

    //
    // Create input channel from samplesheet
    //
    ch_input = ch_samplesheet
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

    //
    // Run simulation subworkflow
    //
    SIMULATE_READS(ch_input)
    ch_versions = ch_versions.mix(SIMULATE_READS.out.versions)

    emit:
    ont_reads        = SIMULATE_READS.out.ont_reads
    pacbio_reads     = SIMULATE_READS.out.pacbio_reads
    illumina_reads   = SIMULATE_READS.out.illumina_reads
    qc_summary_table = SIMULATE_READS.out.qc_summary_table
    versions         = ch_versions
}

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow {

    main:

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    UTILS_NEXTFLOW_PIPELINE (
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // SUBWORKFLOW: Run nf-core/utils_nfcore_pipeline subworkflow
    //
    UTILS_NFCORE_PIPELINE (
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir
    )

    //
    // Read in samplesheet, validate and stage input files
    //
    ch_samplesheet = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)

    //
    // WORKFLOW: Run main workflow
    //
    SEQAUDIT (
        ch_samplesheet
    )

    emit:
    ont_reads        = SEQAUDIT.out.ont_reads
    pacbio_reads     = SEQAUDIT.out.pacbio_reads
    illumina_reads   = SEQAUDIT.out.illumina_reads
    qc_summary_table = SEQAUDIT.out.qc_summary_table
}
