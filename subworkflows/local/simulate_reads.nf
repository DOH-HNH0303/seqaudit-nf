include { FETCH_GENOME } from '../../modules/local/fetch_genome'
include { PBSIM3_ONT_MULTI } from '../../modules/local/pbsim3_ont_multi'
include { PBSIM3_PACBIO_MULTI } from '../../modules/local/pbsim3_pacbio_multi'
include { ART_ILLUMINA_MULTI } from '../../modules/local/art_illumina_multi'
include { FASTQ_QC_CONSOLIDATED } from '../../modules/local/fastq_qc_consolidated'

workflow SIMULATE_READS {
    take:
    ch_input // channel: [meta, genome_id]

    main:
    ch_versions = Channel.empty()
    ch_all_ont_reads = Channel.empty()
    ch_all_pacbio_reads = Channel.empty()
    ch_all_illumina_reads = Channel.empty()

    // Debug: Print parameters
    log.info "ONT simulator: ${params.ont_simulator}"
    log.info "ONT model URL: ${params.ont_model_url}"

    // Fetch genome sequences
    FETCH_GENOME(ch_input)
    ch_versions = ch_versions.mix(FETCH_GENOME.out.versions)

    // Branch based on simulation requirements
    ch_genomes = FETCH_GENOME.out.genome

    // Debug: Check what's in the channel
    ch_genomes.view { meta, genome -> "Genome channel: ${meta}" }

    // ONT simulation branch - generate multiple datasets per sample
    ch_ont_input = ch_genomes.filter { meta, genome ->
        log.info "Checking ONT reads for ${meta.id}: ${meta.ont_reads}"
        meta.ont_reads > 0
    }

    // Debug: Check filtered channel
    ch_ont_input.view { meta, genome -> "ONT input: ${meta.id} with ${meta.ont_reads} datasets" }

    if (params.ont_simulator == 'pbsim3') {
        log.info "Using PBSIM3 for ONT simulation"
        ont_model_file = Channel.fromPath(params.ont_model_url)
        PBSIM3_ONT_MULTI(ch_ont_input, ont_model_file)
        ch_all_ont_reads = PBSIM3_ONT_MULTI.out.reads
            .map { meta, reads -> reads }
            .flatten()
            .collect()
        ch_versions = ch_versions.mix(PBSIM3_ONT_MULTI.out.versions)
    } else {
        log.info "No ONT simulator specified or unrecognized: ${params.ont_simulator}"
        ch_all_ont_reads = Channel.empty()
    }

    // PacBio simulation branch - generate multiple datasets per sample
    ch_pacbio_input = ch_genomes.filter { meta, genome -> meta.pacbio_reads > 0 }

    if (params.pacbio_simulator == 'pbsim3') {
        pacbio_model_file = Channel.fromPath(params.pacbio_model_url)
        PBSIM3_PACBIO_MULTI(ch_pacbio_input, pacbio_model_file)
        ch_all_pacbio_reads = PBSIM3_PACBIO_MULTI.out.reads
            .map { meta, reads -> reads }
            .flatten()
            .collect()
        ch_versions = ch_versions.mix(PBSIM3_PACBIO_MULTI.out.versions)
    } else {
        ch_all_pacbio_reads = Channel.empty()
    }

    // Illumina simulation branch - generate multiple datasets per sample
    ch_illumina_input = ch_genomes.filter { meta, genome -> meta.illumina_reads > 0 }

    ART_ILLUMINA_MULTI(ch_illumina_input)
    ch_all_illumina_reads = ART_ILLUMINA_MULTI.out.reads
        .map { meta, r1_files, r2_files -> [r1_files, r2_files] }
        .transpose()
        .flatten()
        .collect()
    ch_versions = ch_versions.mix(ART_ILLUMINA_MULTI.out.versions)

    // Create consolidated QC report
    ch_versions_yml = ch_versions
        .unique()
        .collectFile(name: 'software_versions.yml')

    FASTQ_QC_CONSOLIDATED(
        ch_all_ont_reads.ifEmpty([]),
        ch_all_pacbio_reads.ifEmpty([]),
        ch_all_illumina_reads.ifEmpty([]),
        ch_versions_yml
    )

    emit:
    ont_reads = ch_all_ont_reads
    pacbio_reads = ch_all_pacbio_reads
    illumina_reads = ch_all_illumina_reads
    ont_qc_stats = FASTQ_QC_CONSOLIDATED.out.ont_stats
    pacbio_qc_stats = FASTQ_QC_CONSOLIDATED.out.pacbio_stats
    illumina_qc_stats = FASTQ_QC_CONSOLIDATED.out.illumina_stats
    versions = ch_versions
}
