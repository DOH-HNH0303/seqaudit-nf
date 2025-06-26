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

    // Debug: Print parameters with timestamp
    def startTime = System.currentTimeMillis()
    log.info "🚀 [${new Date()}] SIMULATE_READS workflow started"
    log.info "📋 ONT simulator: ${params.ont_simulator}"
    log.info "📋 ONT model URL: ${params.ont_model_url}"
    log.info "📋 PacBio simulator: ${params.pacbio_simulator}"
    log.info "📋 PacBio model URL: ${params.pacbio_model_url}"

    // Debug: Log input channel contents
    ch_input.view { meta, genome_id ->
        "🔍 [${System.currentTimeMillis() - startTime}ms] INPUT: ${meta.id} - ont:${meta.ont_reads}, pacbio:${meta.pacbio_reads}, illumina:${meta.illumina_reads}, source:${meta.genome_source}, id:${genome_id}"
    }

    // Fetch genome sequences
    log.info "🧬 [${System.currentTimeMillis() - startTime}ms] Starting FETCH_GENOME process"
    FETCH_GENOME(ch_input)
    ch_versions = ch_versions.mix(FETCH_GENOME.out.versions)

    // Branch based on simulation requirements
    ch_genomes = FETCH_GENOME.out.genome

    // Debug: Check what's in the genome channel with timing
    ch_genomes.view { meta, genome ->
        "🧬 [${System.currentTimeMillis() - startTime}ms] GENOME_FETCHED: ${meta.id} - ont:${meta.ont_reads}, pacbio:${meta.pacbio_reads}, illumina:${meta.illumina_reads}, genome_size:${genome.size()} bytes"
    }

    // Create separate channels for each technology using direct filtering
    log.info "🔀 [${System.currentTimeMillis() - startTime}ms] Creating filtered channels for each technology"
    //ch_genomes.view()
    ch_input.map { meta, genome ->tuple(*([
        meta.id,
        meta.genome_source,
        meta.genome_id,
        meta.ont_reads,
        meta.pacbio_reads,
        meta.illumina_reads,
        genome]))
    }.collect(flat: false).set
  //.view()

    // ch_genomes.map{meta, genome -> [meta], genome}.spread().collect().view()
    ch_ont_input = ch_genomes.filter { meta, genome ->
        def willProcess = meta.ont_reads > 0
        log.info "🔍 [${System.currentTimeMillis() - startTime}ms] ONT_FILTER: ${meta.id} - reads:${meta.ont_reads}, will_process:${willProcess}"
        return willProcess
    }

    ch_pacbio_input = ch_genomes.filter { meta, genome ->
        def willProcess = meta.pacbio_reads > 0
        log.info "🔍 [${System.currentTimeMillis() - startTime}ms] PACBIO_FILTER: ${meta.id} - reads:${meta.pacbio_reads}, will_process:${willProcess}"
        return willProcess
    }
    // Debug: Check filtered results
    ch_ont_input.view { meta, genome -> "DEBUG: ONT filtered: ${meta.id}" }
    ch_pacbio_input.view { meta, genome -> "DEBUG: PacBio filtered: ${meta.id}" }

    ch_illumina_input = ch_genomes.filter { meta, genome ->
        def willProcess = meta.illumina_reads > 0
        log.info "🔍 [${System.currentTimeMillis() - startTime}ms] ILLUMINA_FILTER: ${meta.id} - reads:${meta.illumina_reads}, will_process:${willProcess}"
        return willProcess
    }

    // Debug: Check filtered channels with detailed timing
    ch_ont_input.view { meta, genome ->
        "✅ [${System.currentTimeMillis() - startTime}ms] ONT_FILTERED_INPUT: ${meta.id} with ${meta.ont_reads} datasets, genome:${genome.name}"
    }

    ch_pacbio_input.view { meta, genome ->
        "✅ [${System.currentTimeMillis() - startTime}ms] PACBIO_FILTERED_INPUT: ${meta.id} with ${meta.pacbio_reads} datasets, genome:${genome.name}"
    }

    ch_illumina_input.view { meta, genome ->
        "✅ [${System.currentTimeMillis() - startTime}ms] ILLUMINA_FILTERED_INPUT: ${meta.id} with ${meta.illumina_reads} datasets, genome:${genome.name}"
    }

    // ONT simulation branch - generate multiple datasets per sample
    log.info "🧪 [${System.currentTimeMillis() - startTime}ms] Starting ONT simulation branch"

    if (params.ont_simulator == 'pbsim3') {
        log.info "🔬 [${System.currentTimeMillis() - startTime}ms] Using PBSIM3 for ONT simulation"
        ont_model_file = Channel.fromPath(params.ont_model_url)

        // Debug: Count items in ONT input channel before processing
        ch_ont_input.count().view { count ->
            "📊 [${System.currentTimeMillis() - startTime}ms] ONT_INPUT_COUNT: ${count} samples will be processed"
        }

        PBSIM3_ONT_MULTI(ch_ont_input, ont_model_file)

        // FIXED: Proper channel handling for ONT reads with debug
        ch_all_ont_reads = PBSIM3_ONT_MULTI.out.reads
            .map { meta, reads ->
                log.info "🧬 [${System.currentTimeMillis() - startTime}ms] ONT_OUTPUT: ${meta.id} produced ${reads instanceof List ? reads.size() : 1} read files"
                // Ensure reads is a list and flatten if needed
                def readsList = reads instanceof List ? reads.flatten() : [reads]
                return [meta, readsList]
            }
        ch_versions = ch_versions.mix(PBSIM3_ONT_MULTI.out.versions)
    } else {
        log.info "⚠️ [${System.currentTimeMillis() - startTime}ms] No ONT simulator specified or unrecognized: ${params.ont_simulator}"
        ch_all_ont_reads = Channel.empty()
    }

    // PacBio simulation branch - generate multiple datasets per sample
    log.info "🧪 [${System.currentTimeMillis() - startTime}ms] Starting PacBio simulation branch"
    log.info "🔬 [${System.currentTimeMillis() - startTime}ms] Using PBSIM3 for PacBio simulation"
    pacbio_model_file = Channel.fromPath(params.pacbio_model_url)

    // Debug: Count items in PacBio input channel before processing
    ch_pacbio_input.count().view { count ->
        "📊 [${System.currentTimeMillis() - startTime}ms] PACBIO_INPUT_COUNT: ${count} samples will be processed"
    }

    PBSIM3_PACBIO_MULTI(ch_pacbio_input, pacbio_model_file)

    // FIXED: Proper channel handling for PacBio reads with debug
    ch_all_pacbio_reads = PBSIM3_PACBIO_MULTI.out.reads
        .map { meta, reads ->
            log.info "🧬 [${System.currentTimeMillis() - startTime}ms] PACBIO_OUTPUT: ${meta.id} produced ${reads instanceof List ? reads.size() : 1} read files"
            // Ensure reads is a list and flatten if needed
            def readsList = reads instanceof List ? reads.flatten() : [reads]
            return [meta, readsList]
        }
    ch_versions = ch_versions.mix(PBSIM3_PACBIO_MULTI.out.versions)

    // if (params.pacbio_simulator == 'pbsim3') {
    //     log.info "🔬 [${System.currentTimeMillis() - startTime}ms] Using PBSIM3 for PacBio simulation"
    //     pacbio_model_file = Channel.fromPath(params.pacbio_model_url)

    //     // Debug: Count items in PacBio input channel before processing
    //     ch_pacbio_input.count().view { count ->
    //         "📊 [${System.currentTimeMillis() - startTime}ms] PACBIO_INPUT_COUNT: ${count} samples will be processed"
    //     }

    //     PBSIM3_PACBIO_MULTI(ch_pacbio_input, pacbio_model_file)

    //     // FIXED: Proper channel handling for PacBio reads with debug
    //     ch_all_pacbio_reads = PBSIM3_PACBIO_MULTI.out.reads
    //         .map { meta, reads ->
    //             log.info "🧬 [${System.currentTimeMillis() - startTime}ms] PACBIO_OUTPUT: ${meta.id} produced ${reads instanceof List ? reads.size() : 1} read files"
    //             // Ensure reads is a list and flatten if needed
    //             def readsList = reads instanceof List ? reads.flatten() : [reads]
    //             return [meta, readsList]
    //         }
    //     ch_versions = ch_versions.mix(PBSIM3_PACBIO_MULTI.out.versions)
    // } else {
    //     log.info "⚠️ [${System.currentTimeMillis() - startTime}ms] No PacBio simulator specified or unrecognized: ${params.pacbio_simulator}"
    //     ch_all_pacbio_reads = Channel.empty()
    // }

    // Illumina simulation branch - generate multiple datasets per sample
    log.info "🧪 [${System.currentTimeMillis() - startTime}ms] Starting Illumina simulation branch"

    // Debug: Count items in Illumina input channel before processing
    ch_illumina_input.count().view { count ->
        "📊 [${System.currentTimeMillis() - startTime}ms] ILLUMINA_INPUT_COUNT: ${count} samples will be processed"
    }

    // Process Illumina samples
    ART_ILLUMINA_MULTI(ch_illumina_input)
    ch_all_illumina_reads = ART_ILLUMINA_MULTI.out.reads
        .map { meta, r1_files, r2_files ->
            log.info "🧬 [${System.currentTimeMillis() - startTime}ms] ILLUMINA_OUTPUT: ${meta.id} produced R1:${r1_files instanceof List ? r1_files.size() : 1} R2:${r2_files instanceof List ? r2_files.size() : 1} files"
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

    // Prepare data for manifest creation using a simpler approach
    log.info "📋 [${System.currentTimeMillis() - startTime}ms] Starting manifest creation process"

    // Create maps for each technology type with debug logging
    ch_ont_files = ch_all_ont_reads
        .map { meta, reads ->
            def firstRead = reads instanceof List ? reads[0] : reads
            log.info "📁 [${System.currentTimeMillis() - startTime}ms] ONT_FILE_MAP: ${meta.id} -> ${firstRead?.name ?: 'null'}"
            return [meta.id, firstRead]
        }

    ch_pacbio_files = ch_all_pacbio_reads
        .map { meta, reads ->
            def firstRead = reads instanceof List ? reads[0] : reads
            log.info "📁 [${System.currentTimeMillis() - startTime}ms] PACBIO_FILE_MAP: ${meta.id} -> ${firstRead?.name ?: 'null'}"
            return [meta.id, firstRead]
        }

    ch_illumina_files = ch_all_illumina_reads
        .map { meta, r1_files, r2_files ->
            def firstR1 = r1_files instanceof List ? r1_files[0] : r1_files
            def firstR2 = r2_files instanceof List ? r2_files[0] : r2_files
            log.info "📁 [${System.currentTimeMillis() - startTime}ms] ILLUMINA_FILE_MAP: ${meta.id} -> R1:${firstR1?.name ?: 'null'}, R2:${firstR2?.name ?: 'null'}"
            return [meta.id, firstR1, firstR2]
        }

    // Debug: Count files in each channel
    ch_ont_files.count().view { count ->
        "📊 [${System.currentTimeMillis() - startTime}ms] ONT_FILES_COUNT: ${count}"
    }
    ch_pacbio_files.count().view { count ->
        "📊 [${System.currentTimeMillis() - startTime}ms] PACBIO_FILES_COUNT: ${count}"
    }
    ch_illumina_files.count().view { count ->
        "📊 [${System.currentTimeMillis() - startTime}ms] ILLUMINA_FILES_COUNT: ${count}"
    }

    // Get all sample metadata and create manifest input
    log.info "🔗 [${System.currentTimeMillis() - startTime}ms] Joining channels for manifest creation"

    ch_manifest_input = ch_genomes
        .map { meta, genome ->
            log.info "🔗 [${System.currentTimeMillis() - startTime}ms] GENOME_FOR_MANIFEST: ${meta.id}"
            return [meta.id, meta]
        }
        .join(ch_ont_files, remainder: true)
        .join(ch_pacbio_files, remainder: true)
        .join(ch_illumina_files, remainder: true)
        .map { tuple ->
            def sample_id = tuple[0]
            def meta = tuple[1]
            def ont_file = tuple.size() > 2 && tuple[2] != null ? tuple[2] : null
            def pacbio_file = tuple.size() > 3 && tuple[3] != null ? tuple[3] : null
            def illumina_r1 = tuple.size() > 4 && tuple[4] != null ? tuple[4] : null
            def illumina_r2 = tuple.size() > 5 && tuple[5] != null ? tuple[5] : null

            log.info "🔗 [${System.currentTimeMillis() - startTime}ms] MANIFEST_JOIN: ${sample_id} - ont:${ont_file ? 'present' : 'missing'}, pacbio:${pacbio_file ? 'present' : 'missing'}, illumina:${illumina_r1 && illumina_r2 ? 'present' : 'missing'}"

            // Handle missing files by providing unique placeholder files
            def ont_path = ont_file ? ont_file : file("NO_FILE_ont_${sample_id}")
            def pacbio_path = pacbio_file ? pacbio_file : file("NO_FILE_pacbio_${sample_id}")
            def illumina_r1_path = illumina_r1 ? illumina_r1 : file("NO_FILE_illumina_r1_${sample_id}")
            def illumina_r2_path = illumina_r2 ? illumina_r2 : file("NO_FILE_illumina_r2_${sample_id}")

            return [meta, ont_path, pacbio_path, illumina_r1_path, illumina_r2_path]
        }

    ch_manifest_input.view { meta, ont, pacbio, r1, r2 ->
        "📋 [${System.currentTimeMillis() - startTime}ms] MANIFEST_INPUT: ${meta.id} - ont:${ont.name}, pacbio:${pacbio.name}, r1:${r1.name}, r2:${r2.name}"
    }

    // Create individual manifests for each sample
    log.info "📝 [${System.currentTimeMillis() - startTime}ms] Creating individual manifests"
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
