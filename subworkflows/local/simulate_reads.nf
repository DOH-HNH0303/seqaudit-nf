include { FETCH_GENOME } from '../modules/local/fetch_genome'
include { PBSIM3_ONT } from '../modules/local/pbsim3_ont'
include { PBSIM3_PACBIO } from '../modules/local/pbsim3_pacbio'
include { NANOSIM } from '../modules/local/nanosim'
include { ART_ILLUMINA } from '../modules/local/art_illumina'

workflow SIMULATE_READS {
    take:
    ch_input // channel: [meta, genome_id]

    main:
    ch_versions = Channel.empty()

    // Fetch genome sequences
    FETCH_GENOME(ch_input)
    ch_versions = ch_versions.mix(FETCH_GENOME.out.versions)

    // Branch based on simulation requirements
    ch_genomes = FETCH_GENOME.out.genome

    // ONT simulation branch
    ch_ont_input = ch_genomes.filter { meta, genome -> meta.ont_reads > 0 }
    ch_ont_reads = Channel.empty()

    if (params.ont_simulator == 'pbsim3') {
        PBSIM3_ONT(ch_ont_input)
        ch_ont_reads = PBSIM3_ONT.out.reads
        ch_versions = ch_versions.mix(PBSIM3_ONT.out.versions)
    } else if (params.ont_simulator == 'nanosim') {
        NANOSIM(ch_ont_input)
        ch_ont_reads = NANOSIM.out.reads
        ch_versions = ch_versions.mix(NANOSIM.out.versions)
    }

    // PacBio simulation branch
    ch_pacbio_input = ch_genomes.filter { meta, genome -> meta.pacbio_reads > 0 }
    ch_pacbio_reads = Channel.empty()

    if (params.pacbio_simulator == 'pbsim3') {
        PBSIM3_PACBIO(ch_pacbio_input)
        ch_pacbio_reads = PBSIM3_PACBIO.out.reads
        ch_versions = ch_versions.mix(PBSIM3_PACBIO.out.versions)
    }

    // Illumina simulation branch
    ch_illumina_input = ch_genomes.filter { meta, genome -> meta.illumina_reads > 0 }
    ch_illumina_reads = Channel.empty()

    ART_ILLUMINA(ch_illumina_input)
    ch_illumina_reads = ART_ILLUMINA.out.reads
    ch_versions = ch_versions.mix(ART_ILLUMINA.out.versions)

    emit:
    ont_reads = ch_ont_reads
    pacbio_reads = ch_pacbio_reads
    illumina_reads = ch_illumina_reads
    qc_reports = ch_versions
    versions = ch_versions
}
