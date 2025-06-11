process FASTQ_QC_ONT {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::seqkit=2.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.3.1--h9ee0642_0' :
        'quay.io/biocontainers/seqkit:2.3.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}_qc_stats.txt"), emit: stats
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Initialize stats file
    echo "Sample: ${meta.id}" > ${prefix}_qc_stats.txt
    echo "Platform: ONT" >> ${prefix}_qc_stats.txt
    echo "===================" >> ${prefix}_qc_stats.txt
    echo "" >> ${prefix}_qc_stats.txt

    # Process each read file
    for read_file in ${reads}; do
        echo "File: \$(basename \$read_file)" >> ${prefix}_qc_stats.txt
        echo "-------------------" >> ${prefix}_qc_stats.txt

        # Basic statistics using seqkit
        seqkit stats \$read_file >> ${prefix}_qc_stats.txt
        echo "" >> ${prefix}_qc_stats.txt
    done

    # Summary statistics
    echo "Summary:" >> ${prefix}_qc_stats.txt
    echo "--------" >> ${prefix}_qc_stats.txt
    total_reads=\$(seqkit stats ${reads} | tail -n +2 | awk '{sum+=\$4} END {print sum}')
    total_bases=\$(seqkit stats ${reads} | tail -n +2 | awk '{sum+=\$5} END {print sum}')
    avg_length=\$(seqkit stats ${reads} | tail -n +2 | awk '{sum+=\$6*\$4; count+=\$4} END {if(count>0) print sum/count; else print 0}')

    echo "Total reads: \$total_reads" >> ${prefix}_qc_stats.txt
    echo "Total bases: \$total_bases" >> ${prefix}_qc_stats.txt
    echo "Average read length: \$avg_length" >> ${prefix}_qc_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | cut -d' ' -f2)
    END_VERSIONS
    """
}
