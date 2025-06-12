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
    # Remove pipefail temporarily to prevent SIGPIPE issues
    set -euo

    # Function to safely run commands
    safe_run() {
        "\$@" || {
            echo "Command failed: \$*" >&2
            return 1
        }
    }

    # Function to safely count files
    count_files() {
        local pattern="\$1"
        find . -name "\$pattern" -type f 2>/dev/null | wc -l || echo "0"
    }

    # Function to safely list files
    list_files() {
        local pattern="\$1"
        find . -name "\$pattern" -type f 2>/dev/null || true
    }

    # Debug: List all files in working directory
    echo "DEBUG: Files in working directory:"
    ls -la 2>/dev/null || echo "No files found"
    echo "DEBUG: Looking for dataset files:"
    find . -name "*dataset*.fastq.gz" 2>/dev/null || echo "No dataset files found"
    echo "DEBUG: End of file listing"

    # Initialize consolidated QC files
    echo "=== CONSOLIDATED QC STATISTICS ===" > consolidated_qc_summary.txt
    echo "Generated on: \$(date)" >> consolidated_qc_summary.txt
    echo "" >> consolidated_qc_summary.txt

    # Process ONT reads if they exist
    ont_files=\$(count_files "*ont_dataset*.fastq.gz")
    echo "DEBUG: Found \$ont_files ONT files"

    if [ "\$ont_files" -gt 0 ]; then
        echo "ONT SEQUENCING DATA QUALITY CONTROL" > ont_qc_stats.txt
        echo "====================================" >> ont_qc_stats.txt
        echo "Generated on: \$(date)" >> ont_qc_stats.txt
        echo "" >> ont_qc_stats.txt

        echo "Individual Dataset Statistics:" >> ont_qc_stats.txt
        echo "-----------------------------" >> ont_qc_stats.txt

        # Process each file individually
        list_files "*ont_dataset*.fastq.gz" | while IFS= read -r read_file; do
            if [ -n "\$read_file" ] && [ -f "\$read_file" ]; then
                echo "" >> ont_qc_stats.txt
                echo "File: \$read_file" >> ont_qc_stats.txt
                echo "-------------------" >> ont_qc_stats.txt
                if command -v seqkit >/dev/null 2>&1; then
                    seqkit stats "\$read_file" >> ont_qc_stats.txt 2>/dev/null || echo "Error processing \$read_file" >> ont_qc_stats.txt
                else
                    echo "seqkit not available" >> ont_qc_stats.txt
                fi
            fi
        done

        echo "" >> ont_qc_stats.txt
        echo "CONSOLIDATED ONT SUMMARY:" >> ont_qc_stats.txt
        echo "========================" >> ont_qc_stats.txt

        # Get list of files
        ont_file_list=\$(list_files "*ont_dataset*.fastq.gz")
        if [ -n "\$ont_file_list" ]; then
            # Process files with seqkit if available
            if command -v seqkit >/dev/null 2>&1; then
                echo "\$ont_file_list" | tr '\\n' '\\0' | xargs -0 seqkit stats >> ont_qc_stats.txt 2>/dev/null || echo "Error in consolidated ONT stats" >> ont_qc_stats.txt

                # Calculate totals with error handling
                total_ont_reads=\$(echo "\$ont_file_list" | tr '\\n' '\\0' | xargs -0 seqkit stats 2>/dev/null | tail -n +2 | awk '{sum+=\$4} END {print (sum ? sum : 0)}' || echo "0")
                total_ont_bases=\$(echo "\$ont_file_list" | tr '\\n' '\\0' | xargs -0 seqkit stats 2>/dev/null | tail -n +2 | awk '{sum+=\$5} END {print (sum ? sum : 0)}' || echo "0")
                avg_ont_length=\$(echo "\$ont_file_list" | tr '\\n' '\\0' | xargs -0 seqkit stats 2>/dev/null | tail -n +2 | awk '{sum+=\$6*\$4; count+=\$4} END {if(count>0) print sum/count; else print 0}' || echo "0")
            else
                total_ont_reads="N/A"
                total_ont_bases="N/A"
                avg_ont_length="N/A"
            fi

            echo "" >> ont_qc_stats.txt
            echo "TOTAL ONT STATISTICS:" >> ont_qc_stats.txt
            echo "Total reads: \$total_ont_reads" >> ont_qc_stats.txt
            echo "Total bases: \$total_ont_bases" >> ont_qc_stats.txt
            echo "Average read length: \$avg_ont_length" >> ont_qc_stats.txt
            echo "Number of datasets: \$(echo "\$ont_file_list" | wc -l)" >> ont_qc_stats.txt
        else
            echo "No ONT dataset files found" >> ont_qc_stats.txt
        fi
    else
        echo "No ONT files found for processing" > ont_qc_stats.txt
    fi

    # Process PacBio reads if they exist
    pacbio_files=\$(count_files "*pacbio_dataset*.fastq.gz")
    echo "DEBUG: Found \$pacbio_files PacBio files"

    if [ "\$pacbio_files" -gt 0 ]; then
        echo "PACBIO SEQUENCING DATA QUALITY CONTROL" > pacbio_qc_stats.txt
        echo "=======================================" >> pacbio_qc_stats.txt
        echo "Generated on: \$(date)" >> pacbio_qc_stats.txt
        echo "" >> pacbio_qc_stats.txt

        echo "Individual Dataset Statistics:" >> pacbio_qc_stats.txt
        echo "-----------------------------" >> pacbio_qc_stats.txt

        # Process each file individually
        list_files "*pacbio_dataset*.fastq.gz" | while IFS= read -r read_file; do
            if [ -n "\$read_file" ] && [ -f "\$read_file" ]; then
                echo "" >> pacbio_qc_stats.txt
                echo "File: \$read_file" >> pacbio_qc_stats.txt
                echo "-------------------" >> pacbio_qc_stats.txt
                if command -v seqkit >/dev/null 2>&1; then
                    seqkit stats "\$read_file" >> pacbio_qc_stats.txt 2>/dev/null || echo "Error processing \$read_file" >> pacbio_qc_stats.txt
                else
                    echo "seqkit not available" >> pacbio_qc_stats.txt
                fi
            fi
        done

        echo "" >> pacbio_qc_stats.txt
        echo "CONSOLIDATED PACBIO SUMMARY:" >> pacbio_qc_stats.txt
        echo "============================" >> pacbio_qc_stats.txt

        # Get list of files
        pacbio_file_list=\$(list_files "*pacbio_dataset*.fastq.gz")
        if [ -n "\$pacbio_file_list" ]; then
            # Process files with seqkit if available
            if command -v seqkit >/dev/null 2>&1; then
                echo "\$pacbio_file_list" | tr '\\n' '\\0' | xargs -0 seqkit stats >> pacbio_qc_stats.txt 2>/dev/null || echo "Error in consolidated PacBio stats" >> pacbio_qc_stats.txt

                # Calculate totals with error handling
                total_pacbio_reads=\$(echo "\$pacbio_file_list" | tr '\\n' '\\0' | xargs -0 seqkit stats 2>/dev/null | tail -n +2 | awk '{sum+=\$4} END {print (sum ? sum : 0)}' || echo "0")
                total_pacbio_bases=\$(echo "\$pacbio_file_list" | tr '\\n' '\\0' | xargs -0 seqkit stats 2>/dev/null | tail -n +2 | awk '{sum+=\$5} END {print (sum ? sum : 0)}' || echo "0")
                avg_pacbio_length=\$(echo "\$pacbio_file_list" | tr '\\n' '\\0' | xargs -0 seqkit stats 2>/dev/null | tail -n +2 | awk '{sum+=\$6*\$4; count+=\$4} END {if(count>0) print sum/count; else print 0}' || echo "0")
            else
                total_pacbio_reads="N/A"
                total_pacbio_bases="N/A"
                avg_pacbio_length="N/A"
            fi

            echo "" >> pacbio_qc_stats.txt
            echo "TOTAL PACBIO STATISTICS:" >> pacbio_qc_stats.txt
            echo "Total reads: \$total_pacbio_reads" >> pacbio_qc_stats.txt
            echo "Total bases: \$total_pacbio_bases" >> pacbio_qc_stats.txt
            echo "Average read length: \$avg_pacbio_length" >> pacbio_qc_stats.txt
            echo "Number of datasets: \$(echo "\$pacbio_file_list" | wc -l)" >> pacbio_qc_stats.txt
        else
            echo "No PacBio dataset files found" >> pacbio_qc_stats.txt
        fi
    else
        echo "No PacBio files found for processing" > pacbio_qc_stats.txt
    fi

    # Process Illumina reads if they exist
    illumina_files=\$(count_files "*illumina_dataset*_R1.fastq.gz")
    echo "DEBUG: Found \$illumina_files Illumina R1 files"

    if [ "\$illumina_files" -gt 0 ]; then
        echo "ILLUMINA SEQUENCING DATA QUALITY CONTROL" > illumina_qc_stats.txt
        echo "=========================================" >> illumina_qc_stats.txt
        echo "Generated on: \$(date)" >> illumina_qc_stats.txt
        echo "" >> illumina_qc_stats.txt

        echo "Individual Dataset Statistics:" >> illumina_qc_stats.txt
        echo "-----------------------------" >> illumina_qc_stats.txt

        # Process each file individually
        list_files "*illumina_dataset*_R1.fastq.gz" | while IFS= read -r r1_file; do
            if [ -n "\$r1_file" ] && [ -f "\$r1_file" ]; then
                r2_file=\$(echo "\$r1_file" | sed 's/_R1.fastq.gz/_R2.fastq.gz/')
                dataset_name=\$(echo "\$r1_file" | sed 's/_R1.fastq.gz//')

                echo "" >> illumina_qc_stats.txt
                echo "Dataset: \$dataset_name" >> illumina_qc_stats.txt
                echo "-------------------" >> illumina_qc_stats.txt
                echo "R1 File: \$r1_file" >> illumina_qc_stats.txt

                if command -v seqkit >/dev/null 2>&1; then
                    seqkit stats "\$r1_file" >> illumina_qc_stats.txt 2>/dev/null || echo "Error processing \$r1_file" >> illumina_qc_stats.txt
                else
                    echo "seqkit not available" >> illumina_qc_stats.txt
                fi

                if [ -f "\$r2_file" ]; then
                    echo "" >> illumina_qc_stats.txt
                    echo "R2 File: \$r2_file" >> illumina_qc_stats.txt
                    if command -v seqkit >/dev/null 2>&1; then
                        seqkit stats "\$r2_file" >> illumina_qc_stats.txt 2>/dev/null || echo "Error processing \$r2_file" >> illumina_qc_stats.txt
                    else
                        echo "seqkit not available" >> illumina_qc_stats.txt
                    fi
                else
                    echo "" >> illumina_qc_stats.txt
                    echo "R2 File: \$r2_file (NOT FOUND)" >> illumina_qc_stats.txt
                fi

                # Basic quality information with error handling
                echo "" >> illumina_qc_stats.txt
                echo "Quality Information (first 1000 reads from R1):" >> illumina_qc_stats.txt
                if [ -f "\$r1_file" ] && command -v zcat >/dev/null 2>&1; then
                    qual_lines=\$(zcat "\$r1_file" 2>/dev/null | awk 'NR%4==0' | head -1000 | wc -l 2>/dev/null || echo "0")
                    echo "Quality lines analyzed: \$qual_lines" >> illumina_qc_stats.txt
                    avg_qual_length=\$(zcat "\$r1_file" 2>/dev/null | awk 'NR%4==0' | head -1000 | awk '{sum+=length(\$0); count++} END {if(count>0) print sum/count; else print 0}' 2>/dev/null || echo "0")
                    echo "Average quality string length: \$avg_qual_length" >> illumina_qc_stats.txt
                else
                    echo "Could not analyze quality information - file not accessible or zcat not available" >> illumina_qc_stats.txt
                fi
            fi
        done

        echo "" >> illumina_qc_stats.txt
        echo "CONSOLIDATED ILLUMINA SUMMARY:" >> illumina_qc_stats.txt
        echo "==============================" >> illumina_qc_stats.txt

        # Get list of files
        illumina_file_list=\$(list_files "*illumina_dataset*_R*.fastq.gz")
        if [ -n "\$illumina_file_list" ]; then
            # Process files with seqkit if available
            if command -v seqkit >/dev/null 2>&1; then
                echo "\$illumina_file_list" | tr '\\n' '\\0' | xargs -0 seqkit stats >> illumina_qc_stats.txt 2>/dev/null || echo "Error in consolidated Illumina stats" >> illumina_qc_stats.txt

                # Calculate totals with error handling
                total_illumina_reads=\$(echo "\$illumina_file_list" | tr '\\n' '\\0' | xargs -0 seqkit stats 2>/dev/null | tail -n +2 | awk '{sum+=\$4} END {print (sum ? sum : 0)}' || echo "0")
                total_illumina_bases=\$(echo "\$illumina_file_list" | tr '\\n' '\\0' | xargs -0 seqkit stats 2>/dev/null | tail -n +2 | awk '{sum+=\$5} END {print (sum ? sum : 0)}' || echo "0")
                avg_illumina_length=\$(echo "\$illumina_file_list" | tr '\\n' '\\0' | xargs -0 seqkit stats 2>/dev/null | tail -n +2 | awk '{sum+=\$6*\$4; count+=\$4} END {if(count>0) print sum/count; else print 0}' || echo "0")
            else
                total_illumina_reads="N/A"
                total_illumina_bases="N/A"
                avg_illumina_length="N/A"
            fi

            echo "" >> illumina_qc_stats.txt
            echo "TOTAL ILLUMINA STATISTICS:" >> illumina_qc_stats.txt
            echo "Total reads: \$total_illumina_reads" >> illumina_qc_stats.txt
            echo "Total bases: \$total_illumina_bases" >> illumina_qc_stats.txt
            echo "Average read length: \$avg_illumina_length" >> illumina_qc_stats.txt

            r1_count=\$(count_files "*illumina_dataset*_R1.fastq.gz")
            r_total_count=\$(echo "\$illumina_file_list" | wc -l)
            echo "Number of datasets: \$r1_count" >> illumina_qc_stats.txt
            echo "Number of read files: \$r_total_count" >> illumina_qc_stats.txt
        else
            echo "No Illumina dataset files found" >> illumina_qc_stats.txt
        fi
    else
        echo "No Illumina files found for processing" > illumina_qc_stats.txt
    fi

    # Ensure all output files exist
    touch ont_qc_stats.txt pacbio_qc_stats.txt illumina_qc_stats.txt

    # Create versions file with error handling
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version 2>/dev/null | cut -d' ' -f2 2>/dev/null || echo "unknown")
    END_VERSIONS

    echo "DEBUG: QC process completed successfully"
    """

    stub:
    """
    echo "ONT QC Stats - Stub" > ont_qc_stats.txt
    echo "PacBio QC Stats - Stub" > pacbio_qc_stats.txt
    echo "Illumina QC Stats - Stub" > illumina_qc_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version 2>/dev/null | cut -d' ' -f2 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
