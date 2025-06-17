process FETCH_GENOME {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::entrez-direct=16.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/entrez-direct:16.2--he881be0_1' :
        'quay.io/biocontainers/entrez-direct:24.0--he881be0_0' }"

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
        echo "Processing ${meta.genome_source} genome: ${genome_id}"

        # Define output file based on database selection
        if [[ "${meta.genome_source}" == "refseq" ]]; then
            filter="^>CP|^>NC|^>NZ"  # RefSeq identifiers (including NZ for RefSeq)
        elif [[ "${meta.genome_source}" == "genbank" ]]; then
            filter="^>NZ|^>GCA|^>CP|^>NC"  # GenBank identifiers (more inclusive)
        else
            echo "Invalid database selection: choose 'refseq' or 'genbank'"
            exit 1
        fi

        # Run esearch, elink, efetch with filtering
        echo "Running esearch for assembly: ${genome_id}"
        esearch -db assembly -query "${genome_id}" | \
            elink -target nuccore | \
            efetch -format fasta > "${prefix}_temp.fasta"

        # Check if we got any sequences
        if [[ ! -s "${prefix}_temp.fasta" ]]; then
            echo "ERROR: No sequences downloaded for ${genome_id}"
            echo "This might be due to:"
            echo "1. Invalid accession number"
            echo "2. Network issues"
            echo "3. NCBI database issues"
            exit 1
        fi

        # Apply filtering
        awk -v pat="\$filter" '/^>/ {flag=0} \$0 ~ pat {print; flag=1; next} flag {print}' "${prefix}_temp.fasta" > "${prefix}.fasta"

        # If filtered file is empty, use the original (less strict filtering)
        if [[ ! -s "${prefix}.fasta" ]]; then
            echo "Warning: No sequences matched the filter pattern. Using all downloaded sequences."
            cp "${prefix}_temp.fasta" "${prefix}.fasta"
        fi

        # Clean up temp file
        rm -f "${prefix}_temp.fasta"

        echo "Downloaded ${meta.genome_source} sequences into ${prefix}.fasta"
        echo "Final file size: \$(wc -l < ${prefix}.fasta) lines"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            entrez-direct: \$(esearch -version 2>&1 | head -n1 | sed 's/^.*version //')
        END_VERSIONS
        """
    } else {
        error "Unknown genome source: ${meta.genome_source}"
    }
}
