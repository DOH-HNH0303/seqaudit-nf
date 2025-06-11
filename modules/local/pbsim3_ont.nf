process PBSIM3_ONT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::pbsim3=3.0.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pbsim3:3.0.5--h43eeafb_0' :
        'quay.io/biocontainers/pbsim3:3.0.5--h9948957_0' }"

    input:
    tuple val(meta), path(genome)
    path(model_file)  // Direct model file input

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    path "*.maf.gz", emit: alignments, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def num_reads = meta.ont_reads
    def mean_length = params.ont_read_length_mean
    def sd_length = params.ont_read_length_sd
    def accuracy = params.ont_accuracy

    """
    # Verify model file exists
    if [ ! -f "${model_file}" ]; then
        echo "ERROR: Model file ${model_file} not found!"
        exit 1
    fi

    echo "Using ONT model file: ${model_file}"

    pbsim \\
        --strategy wgs \\
        --method qshmm \\
        --qshmm ${model_file} \\
        --depth ${num_reads} \\
        --genome ${genome} \\
        --prefix ${prefix}_ont \\
        --length-mean ${mean_length} \\
        --length-sd ${sd_length} \\
        --accuracy-mean ${accuracy} \\
        ${args}

    # Compress output files
    gzip *.fastq
    if ls *.maf 1> /dev/null 2>&1; then
        gzip *.maf
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbsim3: \$(pbsim --version 2>&1 | grep -o 'PBSIM [0-9.]*' | sed 's/PBSIM //')
    END_VERSIONS
    """
}
