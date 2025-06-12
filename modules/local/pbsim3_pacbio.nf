process PBSIM3_PACBIO {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::pbsim3=3.0.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pbsim3:3.0.5--h9948957_0' :
        'quay.io/biocontainers/pbsim3:3.0.5--h9948957_0' }"

    input:
    tuple val(meta), path(genome)
    path model_file

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def pacbio_cov_depth = params.pacbio_cov_depth ?: 30
    def read_length_mean = params.pacbio_read_length_mean ?: 8000
    def read_length_sd = params.pacbio_read_length_sd ?: 3000
    def accuracy = params.pacbio_accuracy ?: 0.87

    """
    # Verify model file exists
    if [ ! -f "${model_file}" ]; then
        echo "ERROR: Model file ${model_file} not found!"
        exit 1
    fi

    echo "Using PacBio model file: ${model_file}"

    pbsim \\
        --strategy wgs \\
        --method qshmm \\
        --qshmm ${model_file} \\
        --depth ${pacbio_cov_depth} \\
        --genome ${genome} \\
        --prefix ${prefix}_pacbio \\
        --length-mean ${read_length_mean} \\
        --length-sd ${read_length_sd} \\
        --accuracy-mean ${accuracy} \\
        ${args}

    # List all generated files for debugging
    echo "Generated files:"
    ls -la

    # PBSIM3 generates files like: prefix_0001.fq.gz, prefix_0002.fq.gz, etc.
    # Check for .fq.gz files first (compressed output)
    if ls ${prefix}_pacbio_*.fq.gz 1> /dev/null 2>&1; then
        echo "Found compressed PBSIM3 output files (.fq.gz)"
        # Rename .fq.gz to .fastq.gz for consistency
        for file in ${prefix}_pacbio_*.fq.gz; do
            newname=\$(echo "\$file" | sed 's/\\.fq\\.gz/.fastq.gz/')
            mv "\$file" "\$newname"
            echo "Renamed \$file to \$newname"
        done
    # Check for uncompressed .fq files
    elif ls ${prefix}_pacbio_*.fq 1> /dev/null 2>&1; then
        echo "Found uncompressed PBSIM3 output files (.fq)"
        # Rename and compress .fq to .fastq.gz
        for file in ${prefix}_pacbio_*.fq; do
            newname=\$(echo "\$file" | sed 's/\\.fq/.fastq/')
            mv "\$file" "\$newname"
            gzip "\$newname"
            echo "Renamed and compressed \$file to \$newname.gz"
        done
    else
        echo "ERROR: No FASTQ files found!"
        echo "Available files:"
        ls -la
        echo "Expected pattern: ${prefix}_pacbio_*.fq or ${prefix}_pacbio_*.fq.gz"
        exit 1
    fi

    # Verify we have output files
    if ! ls *.fastq.gz 1> /dev/null 2>&1; then
        echo "ERROR: No .fastq.gz files found after processing!"
        ls -la
        exit 1
    fi

    echo "Final output files:"
    ls -la *.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbsim3: \$(pbsim --version 2>&1 | grep -o 'PBSIM [0-9.]*' | sed 's/PBSIM //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_pacbio_0001.fastq.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbsim3: \$(pbsim --version 2>&1 | grep -o 'PBSIM [0-9.]*' | sed 's/PBSIM //')
    END_VERSIONS
    """
}
