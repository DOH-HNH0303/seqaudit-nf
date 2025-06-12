process FASTQ_QC_CONSOLIDATED {
    label 'process_low'

    conda "bioconda::seqkit=2.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.3.1--h9ee0642_0' :
        'quay.io/biocontainers/seqkit:2.3.1--h9ee0642_0' }"

    publishDir "${params.outdir}/qc_consolidated", mode: 'copy'

    input:
    path ont_reads
    path pacbio_reads  
    path illumina_reads
    path "software_versions.yml"

    output:
    path "ont_qc_stats.txt", emit: ont_stats, optional: true
    path "pacbio_qc_stats.txt", emit: pacbio_stats, optional: true
    path "illumina_qc_stats.txt", emit: illumina_stats, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    # Initialize consolidated QC files
    echo "=== CONSOLIDATED QC STATISTICS ===" > consolidated_qc_summary.txt
    echo "Generated on: \$(date)" >> consolidated_qc_summary.txt
    echo "" >> consolidated_qc_summary.txt

    # Process ONT reads if they exist
    ont_files=\$(find . -name "*ont_dataset*.fastq.gz" 2>/dev/null | wc -l)
    if [ \$ont_files -gt 0 ]; then
        echo "ONT SEQUENCING DATA QUALITY CONTROL" > ont_qc_stats.txt
        echo "====================================" >> ont_qc_stats.txt
        echo "Generated on: \$(date)" >> ont_qc_stats.txt
        echo "" >> ont_qc_stats.txt
        
        echo "Individual Dataset Statistics:" >> ont_qc_stats.txt
        echo "-----------------------------" >> ont_qc_stats.txt
        
        for read_file in *ont_dataset*.fastq.gz; do
            if [ -f "\$read_file" ]; then
                echo "" >> ont_qc_stats.txt
                echo "File: \$read_file" >> ont_qc_stats.txt
                echo "-------------------" >> ont_qc_stats.txt
                seqkit stats "\$read_file" >> ont_qc_stats.txt
            fi
        done
        
        echo "" >> ont_qc_stats.txt
        echo "CONSOLIDATED ONT SUMMARY:" >> ont_qc_stats.txt
        echo "========================" >> ont_qc_stats.txt
        seqkit stats *ont_dataset*.fastq.gz >> ont_qc_stats.txt
        
        # Calculate totals
        total_ont_reads=\$(seqkit stats *ont_dataset*.fastq.gz | tail -n +2 | awk '{sum+=\$4} END {print sum}')
        total_ont_bases=\$(seqkit stats *ont_dataset*.fastq.gz | tail -n +2 | awk '{sum+=\$5} END {print sum}')
        avg_ont_length=\$(seqkit stats *ont_dataset*.fastq.gz | tail -n +2 | awk '{sum+=\$6*\$4; count+=\$4} END {if(count>0) print sum/count; else print 0}')
        
        echo "" >> ont_qc_stats.txt
        echo "TOTAL ONT STATISTICS:" >> ont_qc_stats.txt
        echo "Total reads: \$total_ont_reads" >> ont_qc_stats.txt
        echo "Total bases: \$total_ont_bases" >> ont_qc_stats.txt
        echo "Average read length: \$avg_ont_length" >> ont_qc_stats.txt
        echo "Number of datasets: \$(ls *ont_dataset*.fastq.gz | wc -l)" >> ont_qc_stats.txt
    else
        echo "No ONT files found for processing" > ont_qc_stats.txt
    fi

    # Process PacBio reads if they exist
    pacbio_files=\$(find . -name "*pacbio_dataset*.fastq.gz" 2>/dev/null | wc -l)
    if [ \$pacbio_files -gt 0 ]; then
        echo "PACBIO SEQUENCING DATA QUALITY CONTROL" > pacbio_qc_stats.txt
        echo "=======================================" >> pacbio_qc_stats.txt
        echo "Generated on: \$(date)" >> pacbio_qc_stats.txt
        echo "" >> pacbio_qc_stats.txt
        
        echo "Individual Dataset Statistics:" >> pacbio_qc_stats.txt
        echo "-----------------------------" >> pacbio_qc_stats.txt
        
        for read_file in *pacbio_dataset*.fastq.gz; do
            if [ -f "\$read_file" ]; then
                echo "" >> pacbio_qc_stats.txt
                echo "File: \$read_file" >> pacbio_qc_stats.txt
                echo "-------------------" >> pacbio_qc_stats.txt
                seqkit stats "\$read_file" >> pacbio_qc_stats.txt
            fi
        done
        
        echo "" >> pacbio_qc_stats.txt
        echo "CONSOLIDATED PACBIO SUMMARY:" >> pacbio_qc_stats.txt
        echo "============================" >> pacbio_qc_stats.txt
        seqkit stats *pacbio_dataset*.fastq.gz >> pacbio_qc_stats.txt
        
        # Calculate totals
        total_pacbio_reads=\$(seqkit stats *pacbio_dataset*.fastq.gz | tail -n +2 | awk '{sum+=\$4} END {print sum}')
        total_pacbio_bases=\$(seqkit stats *pacbio_dataset*.fastq.gz | tail -n +2 | awk '{sum+=\$5} END {print sum}')
        avg_pacbio_length=\$(seqkit stats *pacbio_dataset*.fastq.gz | tail -n +2 | awk '{sum+=\$6*\$4; count+=\$4} END {if(count>0) print sum/count; else print 0}')
        
        echo "" >> pacbio_qc_stats.txt
        echo "TOTAL PACBIO STATISTICS:" >> pacbio_qc_stats.txt
        echo "Total reads: \$total_pacbio_reads" >> pacbio_qc_stats.txt
        echo "Total bases: \$total_pacbio_bases" >> pacbio_qc_stats.txt
        echo "Average read length: \$avg_pacbio_length" >> pacbio_qc_stats.txt
        echo "Number of datasets: \$(ls *pacbio_dataset*.fastq.gz | wc -l)" >> pacbio_qc_stats.txt
    else
        echo "No PacBio files found for processing" > pacbio_qc_stats.txt
    fi

    # Process Illumina reads if they exist
    illumina_files=\$(find . -name "*illumina_dataset*_R1.fastq.gz" 2>/dev/null | wc -l)
    if [ \$illumina_files -gt 0 ]; then
        echo "ILLUMINA SEQUENCING DATA QUALITY CONTROL" > illumina_qc_stats.txt
        echo "=========================================" >> illumina_qc_stats.txt
        echo "Generated on: \$(date)" >> illumina_qc_stats.txt
        echo "" >> illumina_qc_stats.txt
        
        echo "Individual Dataset Statistics:" >> illumina_qc_stats.txt
        echo "-----------------------------" >> illumina_qc_stats.txt
        
        for r1_file in *illumina_dataset*_R1.fastq.gz; do
            if [ -f "\$r1_file" ]; then
                r2_file=\$(echo "\$r1_file" | sed 's/_R1.fastq.gz/_R2.fastq.gz/')
                dataset_name=\$(echo "\$r1_file" | sed 's/_R1.fastq.gz//')
                
                echo "" >> illumina_qc_stats.txt
                echo "Dataset: \$dataset_name" >> illumina_qc_stats.txt
                echo "-------------------" >> illumina_qc_stats.txt
                echo "R1 File: \$r1_file" >> illumina_qc_stats.txt
                seqkit stats "\$r1_file" >> illumina_qc_stats.txt
                echo "" >> illumina_qc_stats.txt
                echo "R2 File: \$r2_file" >> illumina_qc_stats.txt
                seqkit stats "\$r2_file" >> illumina_qc_stats.txt
                
                # Quality score distribution for first 1000 reads
                echo "" >> illumina_qc_stats.txt
                echo "Quality Score Distribution (first 1000 reads from R1):" >> illumina_qc_stats.txt
                zcat "\$r1_file" | awk 'NR%4==0' | head -1000 | awk '{
                    for(i=1;i<=length(\$0);i++) {
                        q=sprintf("%d", ord(substr(\$0,i,1))-33)
                        if(q>=0 && q<=50) qual[q]++
                    }
                }
                function ord(c) {
                    return sprintf("%d", c)
                }
                END {
                    for(i=0;i<=50;i++) if(qual[i]>0) printf "Q%d: %d\\n", i, qual[i]
                }' >> illumina_qc_stats.txt
            fi
        done
        
        echo "" >> illumina_qc_stats.txt
        echo "CONSOLIDATED ILLUMINA SUMMARY:" >> illumina_qc_stats.txt
        echo "==============================" >> illumina_qc_stats.txt
        seqkit stats *illumina_dataset*_R*.fastq.gz >> illumina_qc_stats.txt
        
        # Calculate totals
        total_illumina_reads=\$(seqkit stats *illumina_dataset*_R*.fastq.gz | tail -n +2 | awk '{sum+=\$4} END {print sum}')
        total_illumina_bases=\$(seqkit stats *illumina_dataset*_R*.fastq.gz | tail -n +2 | awk '{sum+=\$5} END {print sum}')
        avg_illumina_length=\$(seqkit stats *illumina_dataset*_R*.fastq.gz | tail -n +2 | awk '{sum+=\$6*\$4; count+=\$4} END {if(count>0) print sum/count; else print 0}')
        
        echo "" >> illumina_qc_stats.txt
        echo "TOTAL ILLUMINA STATISTICS:" >> illumina_qc_stats.txt
        echo "Total reads: \$total_illumina_reads" >> illumina_qc_stats.txt
        echo "Total bases: \$total_illumina_bases" >> illumina_qc_stats.txt
        echo "Average read length: \$avg_illumina_length" >> illumina_qc_stats.txt
        echo "Number of datasets: \$(ls *illumina_dataset*_R1.fastq.gz | wc -l)" >> illumina_qc_stats.txt
        echo "Number of read pairs: \$(ls *illumina_dataset*_R*.fastq.gz | wc -l)" >> illumina_qc_stats.txt
    else
        echo "No Illumina files found for processing" > illumina_qc_stats.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    """
    echo "ONT QC Stats - Stub" > ont_qc_stats.txt
    echo "PacBio QC Stats - Stub" > pacbio_qc_stats.txt  
    echo "Illumina QC Stats - Stub" > illumina_qc_stats.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | cut -d' ' -f2)
    END_VERSIONS
    """
}