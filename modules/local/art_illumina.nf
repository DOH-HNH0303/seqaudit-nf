process ART_ILLUMINA {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::art=2016.06.05"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/art:2016.06.05--h516909a_0' :
        'quay.io/biocontainers/art:3.19.15--1' }"

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

    // Validate and set defaults for all variables
    def num_reads = meta.illumina_reads != null ? meta.illumina_reads as Integer : 1000000
    def read_length = params.illumina_read_length != null ? params.illumina_read_length as Integer : 150
    def fragment_mean = params.illumina_fragment_mean != null ? params.illumina_fragment_mean as Integer : 300
    def fragment_sd = params.illumina_fragment_sd != null ? params.illumina_fragment_sd as Integer : 50
    //def system = params.illumina_system ?: 'HS25'
    def genome_size = 3000000000L
    def illumina_q_score = params.illumina_q_score != null ? params.illumina_q_score as Integer : 42

    // Ensure all values are valid before calculation
    if (num_reads <= 0 || read_length <= 0) {
        error "Invalid values: num_reads=${num_reads}, read_length=${read_length}"
    }

def fold_coverage = (num_reads * read_length * 2) / genome_size

    """
    art_illumina \\
        -p \\
        -i ${genome} \\
        -l ${read_length} \\
        -q ${illumina_q_score} \\
        -f ${fold_coverage} \\
        -m ${fragment_mean} \\
        -s ${fragment_sd} \\
        -o ${prefix}_illumina_ \\
        ${args}

    # Rename output files to standard format
    mv ${prefix}_illumina_1.fq ${prefix}_illumina_R1.fastq
    mv ${prefix}_illumina_2.fq ${prefix}_illumina_R2.fastq

    # Compress output files
    gzip ${prefix}_illumina_R1.fastq
    gzip ${prefix}_illumina_R2.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        art: \$(art_illumina --help 2>&1 | grep -o 'Version [0-9.]*' | sed 's/Version //')
    END_VERSIONS
    """
}
