process FETCH_GENOME {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::entrez-direct=16.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/entrez-direct:16.2--he881be0_1' :
        'quay.io/biocontainers/entrez-direct:16.2--he881be0_1' }"

    input:
    tuple val(meta), val(genome_id)

    output:
    tuple val(meta), path("*.fasta"), emit: genome
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (meta.genome_source == 'local') {
        """
        # Copy local file
        cp ${genome_id} ${prefix}.fasta

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            local_copy: "file copied"
        END_VERSIONS
        """
    } else if (meta.genome_source == 'refseq' || meta.genome_source == 'genbank') {
        """
        # Download from NCBI
        esearch -db assembly -query "${genome_id}" | \\
        elink -target nuccore | \\
        efetch -format fasta > ${prefix}.fasta

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            entrez-direct: \$(esearch -version 2>&1 | head -n1 | sed 's/^.*version //')
        END_VERSIONS
        """
    } else {
        error "Unknown genome source: ${meta.genome_source}"
    }
}
