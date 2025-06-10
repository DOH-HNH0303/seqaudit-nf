process NANOSIM {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::nanosim=3.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanosim:3.1.0--py39h38f01e4_1' :
        'quay.io/biocontainers/nanosim:3.1.0--py39h38f01e4_1' }"

    input:
    tuple val(meta), path(genome)

    output:
    tuple val(meta), path("*.fastq"), emit: reads
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def num_reads = meta.ont_reads

    """
    # Use pre-trained model for ONT simulation
    simulator.py genome \\
        -rg ${genome} \\
        -c \$CONDA_PREFIX/share/nanosim/pre-trained_models/human_NA12878_DNA_FAB49712_guppy/training \\
        -o ${prefix}_nanosim \\
        -n ${num_reads} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanosim: \$(simulator.py --version 2>&1 | grep -o '[0-9.]*')
    END_VERSIONS
    """
}
