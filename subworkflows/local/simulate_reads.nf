include { FETCH_GENOME } from '../../modules/local/fetch_genome'
include { PBSIM3_ONT_MULTI } from '../../modules/local/pbsim3_ont_multi'
include { PBSIM3_PACBIO_MULTI } from '../../modules/local/pbsim3_pacbio_multi'
include { ART_ILLUMINA_MULTI } from '../../modules/local/art_illumina_multi'
include { FASTQ_QC_CONSOLIDATED } from '../../modules/local/fastq_qc_consolidated'
include { CREATE_MANIFEST } from '../../modules/local/create_manifest'

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
    ch_genomes.view { meta, genome -> "Genome channel: ${meta.id} - ont_reads: ${meta.ont_reads}, genome_source: ${meta.genome_source}, genome_id: ${meta.genome_id}" }

    ch_genomes
    .multiMap { meta, genome ->
        ont: meta.ont_reads > 0 ? [meta, genome] : null
        pacbio: meta.pacbio_reads > 0 ? [meta, genome] : null
        illumina: meta.illumina_reads > 0 ? [meta, genome] : null
    }
    .set { ch_branched }

    // Debug: Check branching logic
    ch_branched.ont.view { "ONT branch: ${it}" }
    ch_branched.pacbio.view { "PacBio branch: ${it}" }
    ch_branched.illumina.view { "Illumina branch: ${it}" }

    // Then use the individual branches
    ch_ont_input = ch_branched.ont.filter { it != null }
    ch_pacbio_input = ch_branched.pacbio.filter { it != null }
    ch_illumina_input = ch_branched.illumina.filter { it != null }

    // ONT simulation branch - generate multiple datasets per sample
    // ch_ont_input = ch_genomes.filter { meta, genome ->
    //     log.info "Checking ONT reads for ${meta.id}: ${meta.ont_reads}"
    //     meta.ont_reads > 0
    // }

    // Debug: Check filtered channel
    ch_ont_input.view { meta, genome -> "ONT input: ${meta.id} with ${meta.ont_reads} datasets" }

    if (params.ont_simulator == 'pbsim3') {
        log.info "Using PBSIM3 for ONT simulation"
        ont_model_file = Channel.fromPath(params.ont_model_url)
        PBSIM3_ONT_MULTI(ch_ont_input, ont_model_file)
        // FIXED: Proper channel handling for ONT reads
        ch_all_ont_reads = PBSIM3_ONT_MULTI.out.reads
            .map { meta, reads ->
                // Ensure reads is a list and flatten if needed
                def readsList = reads instanceof List ? reads.flatten() : [reads]
                return [meta, readsList]
            }
        ch_versions = ch_versions.mix(PBSIM3_ONT_MULTI.out.versions)
    } else {
        log.info "No ONT simulator specified or unrecognized: ${params.ont_simulator}"
        ch_all_ont_reads = Channel.empty()
    }

    // PacBio simulation branch - generate multiple datasets per sample
    // ch_pacbio_input = ch_genomes.filter { meta, genome -> meta.pacbio_reads > 0 }
    if (params.pacbio_simulator == 'pbsim3') {
        pacbio_model_file = Channel.fromPath(params.pacbio_model_url)
        PBSIM3_PACBIO_MULTI(ch_pacbio_input, pacbio_model_file)
        // FIXED: Proper channel handling for PacBio reads
        ch_all_pacbio_reads = PBSIM3_PACBIO_MULTI.out.reads
            .map { meta, reads ->
                // Ensure reads is a list and flatten if needed
                def readsList = reads instanceof List ? reads.flatten() : [reads]
                return [meta, readsList]
            }
        ch_versions = ch_versions.mix(PBSIM3_PACBIO_MULTI.out.versions)
    } else {
        ch_all_pacbio_reads = Channel.empty()
    }

    // Illumina simulation branch - generate multiple datasets per sample
    //ch_illumina_input = ch_genomes.filter { meta, genome -> meta.illumina_reads > 0 }

    // Process Illumina samples
    ART_ILLUMINA_MULTI(ch_illumina_input)
    ch_all_illumina_reads = ART_ILLUMINA_MULTI.out.reads
        .map { meta, r1_files, r2_files ->
            // Return separate R1 and R2 channels
            return [meta, r1_files, r2_files]
        }
    ch_versions = ch_versions.mix(ART_ILLUMINA_MULTI.out.versions)

    // Create consolidated QC report
    ch_versions_yml = ch_versions
        .unique()
        .collectFile(name: 'software_versions.yml')

    // Prepare channels for QC (flattened for QC purposes)
    ch_ont_for_qc = ch_all_ont_reads
        .map { meta, reads -> reads instanceof List ? reads.flatten() : [reads] }
        .flatten()
        .collect()
        .ifEmpty([])

    ch_pacbio_for_qc = ch_all_pacbio_reads
        .map { meta, reads -> reads instanceof List ? reads.flatten() : [reads] }
        .flatten()
        .collect()
        .ifEmpty([])

    ch_illumina_for_qc = ch_all_illumina_reads
        .map { meta, r1_files, r2_files ->
            def allFiles = []
            if (r1_files instanceof List) {
                allFiles.addAll(r1_files)
            } else if (r1_files != null) {
                allFiles.add(r1_files)
            }
            if (r2_files instanceof List) {
                allFiles.addAll(r2_files)
            } else if (r2_files != null) {
                allFiles.add(r2_files)
            }
            return allFiles
        }
        .flatten()
        .collect()
        .ifEmpty([])

    FASTQ_QC_CONSOLIDATED(
        ch_ont_for_qc,
        ch_pacbio_for_qc,
        ch_illumina_for_qc,
        ch_versions_yml
    )

    // Prepare data for manifest creation
    // Get all unique sample metadata
    ch_all_samples = ch_genomes
        .map { meta, genome -> meta }
        .unique { it.id }

    // Prepare channels for joining - add sample ID as key and provide default empty files
    ch_ont_for_manifest = ch_all_ont_reads
        .map { meta, reads ->
            def firstRead = reads instanceof List ? reads[0] : reads
            return [meta.id, firstRead]
        }
        .ifEmpty([])

    ch_pacbio_for_manifest = ch_all_pacbio_reads
        .map { meta, reads ->
            def firstRead = reads instanceof List ? reads[0] : reads
            return [meta.id, firstRead]
        }
        .ifEmpty([])

    ch_illumina_for_manifest = ch_all_illumina_reads
        .map { meta, r1_files, r2_files ->
            def firstR1 = r1_files instanceof List ? r1_files[0] : r1_files
            def firstR2 = r2_files instanceof List ? r2_files[0] : r2_files
            return [meta.id, firstR1, firstR2]
        }
        .ifEmpty([])

    // Create manifest input by joining all channels
    ch_manifest_input = ch_all_samples
        .map { meta -> [meta.id, meta] }
        .join(ch_ont_for_manifest, remainder: true)
        .join(ch_pacbio_for_manifest, remainder: true)
        .join(ch_illumina_for_manifest, remainder: true)
        .map { tuple ->
            def sample_id = tuple[0]
            def meta = tuple[1]
            def ont_file = tuple.size() > 2 ? tuple[2] : null
            def pacbio_file = tuple.size() > 3 ? tuple[3] : null
            def illumina_r1 = tuple.size() > 4 ? tuple[4] : null
            def illumina_r2 = tuple.size() > 5 ? tuple[5] : null

            // Handle missing files by providing unique placeholder files
            def ont_path = ont_file ?: file("NO_FILE_ont_${sample_id}")
            def pacbio_path = pacbio_file ?: file("NO_FILE_pacbio_${sample_id}")
            def illumina_r1_path = illumina_r1 ?: file("NO_FILE_illumina_r1_${sample_id}")
            def illumina_r2_path = illumina_r2 ?: file("NO_FILE_illumina_r2_${sample_id}")

            return [meta, ont_path, pacbio_path, illumina_r1_path, illumina_r2_path]
        }
    ch_manifest_input.view()
    // Create individual manifests for each sample
    CREATE_MANIFEST(ch_manifest_input)
    ch_versions = ch_versions.mix(CREATE_MANIFEST.out.versions)

    emit:
    ont_reads = ch_all_ont_reads
    pacbio_reads = ch_all_pacbio_reads
    illumina_reads = ch_all_illumina_reads
    ont_qc_stats = FASTQ_QC_CONSOLIDATED.out.ont_stats
    pacbio_qc_stats = FASTQ_QC_CONSOLIDATED.out.pacbio_stats
    illumina_qc_stats = FASTQ_QC_CONSOLIDATED.out.illumina_stats
    manifest = CREATE_MANIFEST.out.manifest
    versions = ch_versions
}
