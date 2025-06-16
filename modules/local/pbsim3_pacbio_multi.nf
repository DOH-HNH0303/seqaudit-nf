process PBSIM3_PACBIO_MULTI {
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
    tuple val(meta), path("*final*.fastq.gz"), emit: reads
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def num_datasets = meta.pacbio_reads
    def read_length_mean = params.pacbio_read_length_mean ?: 8000
    def read_length_sd = params.pacbio_read_length_sd ?: 3000
    def accuracy = params.pacbio_accuracy ?: 0.87
    def coverage_per_dataset = params.pacbio_cov_depth ?: 30

    """
    # Verify model file exists
    if [ ! -f "${model_file}" ]; then
        echo "ERROR: Model file ${model_file} not found!"
        exit 1
    fi

    echo "Using PacBio model file: ${model_file}"
    echo "Generating ${num_datasets} PacBio datasets for sample ${meta.id}"

    # Generate multiple datasets
    for i in \$(seq 1 ${num_datasets}); do
        dataset_prefix="${prefix}_pacbio_dataset_\${i}"
        echo "Generating PacBio dataset \${i} with prefix \${dataset_prefix}"

        pbsim \\
            --strategy wgs \\
            --method qshmm \\
            --qshmm ${model_file} \\
            --depth ${coverage_per_dataset} \\
            --genome ${genome} \\
            --prefix \${dataset_prefix} \\
            --length-mean ${read_length_mean} \\
            --length-sd ${read_length_sd} \\
            --accuracy-mean ${accuracy} \\
            ${args}
    done

    # List all generated files for debugging
    echo "Generated files:"
    ls -la

    # Process all generated files
    for file in *_pacbio_dataset_*.fq; do
        if [ -f "\$file" ]; then
            newname=\$(echo "\$file" | sed 's/\\.fq/.fastq/')
            mv "\$file" "\$newname"
            gzip "\$newname"
            echo "Processed \$file to \$newname.gz"
        fi
    done

    # Also handle any .fq.gz files
    for file in *_pacbio_dataset_*.fq.gz; do
        if [ -f "\$file" ]; then
            newname=\$(echo "\$file" | sed 's/\\.fq\\.gz/.fastq.gz/')
            mv "\$file" "\$newname"
            echo "Renamed \$file to \$newname"
        fi
    done

    # Verify we have output files
    if ! ls *_pacbio_dataset_*.fastq.gz 1> /dev/null 2>&1; then
        echo "ERROR: No PacBio dataset .fastq.gz files found after processing!"
        ls -la
        exit 1
    fi

    datasets=\$(ls *_pacbio_dataset_*.fastq.gz | awk -F'_' '{print \$4}' | sort -u)
    for dataset in \$datasets; do
        cat ${prefix}_pacbio_dataset_\${dataset}_*.fastq.gz > ${prefix}_pacbio_dataset_final_\${dataset}.fastq.gz; done

    echo "Final pacbio dataset files:"
    ls -la *_pacbio_dataset_final*.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbsim3: \$(pbsim --version 2>&1 | grep -o 'PBSIM [0-9.]*' | sed 's/PBSIM //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def num_datasets = meta.pacbio_reads
    """
    for i in \$(seq 1 ${num_datasets}); do
        touch ${prefix}_pacbio_dataset_\${i}_0001.fastq.gz
    done
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbsim3: \$(pbsim --version 2>&1 | grep -o 'PBSIM [0-9.]*' | sed 's/PBSIM //')
    END_VERSIONS
    """
}
