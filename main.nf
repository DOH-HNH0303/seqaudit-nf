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

// Include modules and subworkflows
include { validateParameters; paramsHelp; paramsSummaryLog } from 'plugin/nf-validation'
include { SIMULATE_READS } from './subworkflows/local/simulate_reads'
include { MULTIQC } from './modules/local/multiqc'

// Print help message if requested
if (params.help) {
    log.info paramsHelp("nextflow run main.nf --input samplesheet.csv --outdir results")
    exit 0
}

// Validate parameters
//validateParameters()

// Print parameter summary
log.info paramsSummaryLog(workflow)

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

    // Collect QC reports for MultiQC
    ch_multiqc_files = SIMULATE_READS.out.qc_reports.collect()

    // Run MultiQC
    MULTIQC(
        ch_multiqc_files,
        [],
        [],
        []
    )

    // Emit outputs
    emit:
    ont_reads = SIMULATE_READS.out.ont_reads
    pacbio_reads = SIMULATE_READS.out.pacbio_reads
    illumina_reads = SIMULATE_READS.out.illumina_reads
    multiqc_report = MULTIQC.out.report
}
