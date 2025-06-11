process PBSIM3_ONT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::pbsim3=3.0.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pbsim3:3.0.5--h43eeafb_0' :
        'quay.io/biocontainers/pbsim3:3.0.5--h9948957_0' }"

    input:
    tuple val(meta), path(genome)

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
    def model = params.ont_model
    def mean_length = params.ont_read_length_mean
    def sd_length = params.ont_read_length_sd
    def accuracy = params.ont_accuracy

    """
    pbsim \\
        --strategy wgs \\
        --method qshmm \\
        --qshmm \$CONDA_PREFIX/share/pbsim3/data/${model} \\
        --depth ${num_reads} \\
        --genome ${genome} \\
        --prefix ${prefix}_ont \\
        --length-mean ${mean_length} \\
        --length-sd ${sd_length} \\
        --accuracy-mean ${accuracy} \\
        ${args}

    # Compress output files
    gzip *.fastq
    if [ -f *.maf ]; then
        gzip *.maf
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbsim3: \$(pbsim --version 2>&1 | grep -o 'PBSIM [0-9.]*' | sed 's/PBSIM //')
    END_VERSIONS
    """
}
