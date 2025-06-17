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
        # Set database type: "refseq" or "genbank"

        # Define output file based on database selection
        if [[ "${meta.genome_source}" == "refseq" ]]; then
            filter="^>CP|^>NC"  # RefSeq identifiers
        elif [[ "${meta.genome_source}" == "genbank" ]]; then
            filter="^>NZ|^>GCA"  # GenBank identifiers
        else
            echo "Invalid database selection: choose 'refseq' or 'genbank'"
            exit 1
        fi

        # Run esearch, elink, efetch with filtering
        esearch -db assembly -query "${genome_id}" | \
            elink -target nuccore | \
            efetch -format fasta | \
            awk -v pat="\$filter" '/^>/ {flag=0} \$0 ~ pat {print; flag=1; next} flag {print}' > "${prefix}.fasta"


        echo "Downloaded ${meta.genome_source} sequences into ${prefix}.fasta"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            entrez-direct: \$(esearch -version 2>&1 | head -n1 | sed 's/^.*version //')
        END_VERSIONS
        """
    } else {
        error "Unknown genome source: ${meta.genome_source}"
    }
}
