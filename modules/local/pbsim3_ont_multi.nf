process PBSIM3_ONT_MULTI {
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
    def num_datasets = meta.ont_reads
    def read_length_mean = params.ont_read_length_mean ?: 10000
    def read_length_sd = params.ont_read_length_sd ?: 5000
    def accuracy = params.ont_accuracy ?: 0.95
    def coverage_per_dataset = params.ont_cov_depth ?: 30

    """
    # Verify model file exists
    if [ ! -f "${model_file}" ]; then
        echo "ERROR: Model file ${model_file} not found!"
        exit 1
    fi

    echo "Using ONT model file: ${model_file}"
    echo "Generating ${num_datasets} ONT datasets for sample ${meta.id}"

    # Generate multiple datasets
    for i in \$(seq 1 ${num_datasets}); do
        dataset_prefix="${prefix}_ont_dataset_\${i}"
        echo "Generating ONT dataset \${i} with prefix \${dataset_prefix}"

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
    for file in *_ont_dataset_*.fq; do
        if [ -f "\$file" ]; then
            newname=\$(echo "\$file" | sed 's/\\.fq/.fastq/')
            mv "\$file" "\$newname"
            gzip "\$newname"
            echo "Processed \$file to \$newname.gz"
        fi
    done

    # Also handle any .fq.gz files
    for file in *_ont_dataset_*.fq.gz; do
        if [ -f "\$file" ]; then
            newname=\$(echo "\$file" | sed 's/\\.fq\\.gz/.fastq.gz/')
            mv "\$file" "\$newname"
            echo "Renamed \$file to \$newname"
        fi
    done

    # Verify we have output files
    if ! ls *_ont_dataset_*.fastq.gz 1> /dev/null 2>&1; then
        echo "ERROR: No ONT dataset .fastq.gz files found after processing!"
        ls -la
        exit 1
    fi

    datasets="$(ls *_ont_dataset_*.fastq.gz | awk -F'_' '{print $4}' | sort -u)"

    for dataset in $datasets; do
        cat ${prefix}_ont_dataset_${dataset}_*.fastq.gz > ${prefix}_ont_dataset_final_${dataset}.fastq.gz; done

    echo "Final ONT dataset files:"
    ls -la *_ont_dataset_final*.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbsim3: \$(pbsim --version 2>&1 | grep -o 'PBSIM [0-9.]*' | sed 's/PBSIM //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def num_datasets = meta.ont_reads
    """
    for i in \$(seq 1 ${num_datasets}); do
        touch ${prefix}_ont_dataset_final_\${i}.fastq.gz
    done
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbsim3: \$(pbsim --version 2>&1 | grep -o 'PBSIM [0-9.]*' | sed 's/PBSIM //')
    END_VERSIONS
    """
}
