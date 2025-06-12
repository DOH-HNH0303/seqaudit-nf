process ART_ILLUMINA_MULTI {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::art=2016.06.05"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/art:2016.06.05--h516909a_0' :
        'quay.io/biocontainers/art:2016.06.05--h869255c_2' }"

    input:
    tuple val(meta), path(genome)

    output:
    tuple val(meta), path("*_R1.fastq.gz"), path("*_R2.fastq.gz"), emit: reads
    path "*.aln", emit: alignments, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def num_datasets = meta.illumina_reads
    def read_length = params.illumina_read_length != null ? params.illumina_read_length as Integer : 150
    def fragment_mean = params.illumina_fragment_mean != null ? params.illumina_fragment_mean as Integer : 300
    def fragment_sd = params.illumina_fragment_sd != null ? params.illumina_fragment_sd as Integer : 50
    def system = params.illumina_system ?: 'MSv3'
    def genome_size = 3000000000L
    def illumina_q_score_min = params.illumina_q_score_min != null ? params.illumina_q_score_min as Integer : 40
    def illumina_q_score_max = params.illumina_q_score_max != null ? params.illumina_q_score_max as Integer : 50
    def reads_per_dataset = params.illumina_reads_per_dataset ?: 1000000
    def fold_coverage = (reads_per_dataset * read_length * 2) / genome_size

    """
    echo "Generating ${num_datasets} Illumina datasets for sample ${meta.id}"

    # Generate multiple datasets
    for i in \$(seq 1 ${num_datasets}); do
        dataset_prefix="${prefix}_illumina_dataset_\${i}"
        echo "Generating Illumina dataset \${i} with prefix \${dataset_prefix}"
        
        art_illumina \\
            -p \\
            -i ${genome} \\
            -l ${read_length} \\
            --minQ ${illumina_q_score_min} \\
            --maxQ ${illumina_q_score_max} \\
            -f ${fold_coverage} \\
            -m ${fragment_mean} \\
            -s ${fragment_sd} \\
            -ss ${system} \\
            -o \${dataset_prefix}_ \\
            ${args}

        # Rename output files to standard format
        mv \${dataset_prefix}_1.fq \${dataset_prefix}_R1.fastq
        mv \${dataset_prefix}_2.fq \${dataset_prefix}_R2.fastq

        # Compress output files
        gzip \${dataset_prefix}_R1.fastq
        gzip \${dataset_prefix}_R2.fastq
    done

    echo "Final Illumina dataset files:"
    ls -la *_illumina_dataset_*_R*.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        art: \$(art_illumina --help 2>&1 | grep -o 'Version [0-9.]*' | sed 's/Version //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def num_datasets = meta.illumina_reads
    """
    for i in \$(seq 1 ${num_datasets}); do
        touch ${prefix}_illumina_dataset_\${i}_R1.fastq.gz
        touch ${prefix}_illumina_dataset_\${i}_R2.fastq.gz
    done
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        art: \$(art_illumina --help 2>&1 | grep -o 'Version [0-9.]*' | sed 's/Version //')
    END_VERSIONS
    """
}