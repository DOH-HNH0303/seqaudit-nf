process ART_ILLUMINA {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::art=2016.06.05"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/art:2016.06.05--h516909a_0' :
        'quay.io/biocontainers/art:2016.06.05--h516909a_0' }"

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
    def num_reads = meta.illumina_reads
    def read_length = params.illumina_read_length
    def fragment_mean = params.illumina_fragment_mean
    def fragment_sd = params.illumina_fragment_sd
    def system = params.illumina_system

    // Calculate fold coverage from number of reads
    def genome_size = 3000000000 // Approximate human genome size, adjust as needed
    def fold_coverage = (num_reads * read_length * 2) / genome_size

    """
    art_illumina \\
        -p \\
        -ss ${system} \\
        -i ${genome} \\
        -l ${read_length} \\
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
