include { FETCH_GENOME } from '../../modules/local/fetch_genome'
include { PBSIM3_ONT } from '../../modules/local/pbsim3_ont'
include { PBSIM3_PACBIO } from '../../modules/local/pbsim3_pacbio'
include { NANOSIM } from '../../modules/local/nanosim'
include { ART_ILLUMINA } from '../../modules/local/art_illumina'

workflow SIMULATE_READS {
    take:
    ch_input // channel: [meta, genome_id]

    main:
    ch_versions = Channel.empty()
    ch_qc_reports = Channel.empty()

    // Fetch genome sequences
    FETCH_GENOME(ch_input)
    ch_versions = ch_versions.mix(FETCH_GENOME.out.versions)

    // Branch based on simulation requirements
    ch_genomes = FETCH_GENOME.out.genome

    // ONT simulation branch
    ch_ont_input = ch_genomes.filter { meta, genome -> meta.ont_reads > 0 }
    ch_ont_reads = Channel.empty()

    if (params.ont_simulator == 'pbsim3') {
        ont_model_file = Channel.fromPath(params.ont_model_url)
        PBSIM3_ONT(ch_ont_input, ont_model_file)
        ch_ont_reads = PBSIM3_ONT.out.reads
        ch_versions = ch_versions.mix(PBSIM3_ONT.out.versions)

        // QC for ONT reads
        FASTQ_QC(ch_ont_reads.map { meta, reads ->
            [meta + [platform: 'ont'], reads]
        })
        ch_qc_reports = ch_qc_reports.mix(FASTQ_QC.out.stats)
        ch_versions = ch_versions.mix(FASTQ_QC.out.versions)

    } else if (params.ont_simulator == 'nanosim') {
        NANOSIM(ch_ont_input)
        ch_ont_reads = NANOSIM.out.reads
        ch_versions = ch_versions.mix(NANOSIM.out.versions)

        // QC for ONT reads
        FASTQ_QC(ch_ont_reads.map { meta, reads ->
            [meta + [platform: 'ont'], reads]
        })
        ch_qc_reports = ch_qc_reports.mix(FASTQ_QC.out.stats)
        ch_versions = ch_versions.mix(FASTQ_QC.out.versions

    }

    // PacBio simulation branch
    ch_pacbio_input = ch_genomes.filter { meta, genome -> meta.pacbio_reads > 0 }
    ch_pacbio_reads = Channel.empty()

    if (params.pacbio_simulator == 'pbsim3') {
        pacbio_model_file = Channel.fromPath(params.pacbio_model_url)
        PBSIM3_PACBIO(ch_pacbio_input, pacbio_model_file)
        ch_pacbio_reads = PBSIM3_PACBIO.out.reads
        ch_versions = ch_versions.mix(PBSIM3_PACBIO.out.versions)

        // QC for PacBio reads
        FASTQ_QC(ch_pacbio_reads.map { meta, reads ->
            [meta + [platform: 'pacbio'], reads]
        })
        ch_qc_reports = ch_qc_reports.mix(FASTQ_QC.out.stats)
        ch_versions = ch_versions.mix(FASTQ_QC.out.versions)

    }

    // Illumina simulation branch
    ch_illumina_input = ch_genomes.filter { meta, genome -> meta.illumina_reads > 0 }
    ch_illumina_reads = Channel.empty()

    ART_ILLUMINA(ch_illumina_input)
    ch_illumina_reads = ART_ILLUMINA.out.reads
    ch_versions = ch_versions.mix(ART_ILLUMINA.out.versions)

    // QC for Illumina reads
    FASTQ_QC(ch_illumina_reads.map { meta, reads ->
        [meta + [platform: 'illumina'], reads]
    })
    ch_qc_reports = ch_qc_reports.mix(FASTQ_QC.out.stats)
    ch_versions = ch_versions.mix(FASTQ_QC.out.versions)

    // Collect all QC reports and create summary table
    ch_all_qc = ch_qc_reports.map { meta, stats -> stats }.collect()
    ch_versions_yml = ch_versions.collectFile(name: 'software_versions.yml')

    QC_SUMMARY(ch_all_qc, ch_versions_yml)

    emit:
    ont_reads = ch_ont_reads
    pacbio_reads = ch_pacbio_reads
    illumina_reads = ch_illumina_reads
    qc_reports = ch_qc_reports
    qc_summary_table = QC_SUMMARY.out.summary_table
    versions = ch_versions
}
