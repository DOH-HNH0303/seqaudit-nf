process QC_SUMMARY {
    label 'process_single'

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    publishDir "${params.outdir}/qc_summary", mode: 'copy'

    input:
    path qc_files
    path "software_versions.yml"

    output:
    path "qc_summary_table.tsv", emit: summary_table
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    python3 << 'EOF'
import os
import re
from pathlib import Path

# Initialize data structure
qc_data = {}

# Process each QC stats file
qc_files = [f for f in os.listdir('.') if f.endswith('_qc_stats.txt')]

for qc_file in qc_files:
    sample_name = qc_file.replace('_qc_stats.txt', '')
    qc_data[sample_name] = {
        'sample': sample_name,
        'platform': 'unknown',
        'total_reads': 0,
        'total_bases': 0,
        'avg_length': 0,
        'files': [],
        'q30_plus': 0,
        'q40_plus': 0
    }

    with open(qc_file, 'r') as f:
        content = f.read()

        # Extract basic stats
        total_reads_match = re.search(r'Total reads: (\\d+)', content)
        if total_reads_match:
            qc_data[sample_name]['total_reads'] = int(total_reads_match.group(1))

        total_bases_match = re.search(r'Total bases: ([\\d,]+)', content)
        if total_bases_match:
            qc_data[sample_name]['total_bases'] = int(total_bases_match.group(1).replace(',', ''))

        avg_length_match = re.search(r'Average read length: ([\\d.]+)', content)
        if avg_length_match:
            qc_data[sample_name]['avg_length'] = float(avg_length_match.group(1))

        # Extract file information
        file_matches = re.findall(r'File: ([^\\n]+)', content)
        qc_data[sample_name]['files'] = file_matches

        # Determine platform based on file names
        if any('illumina' in f.lower() for f in file_matches):
            qc_data[sample_name]['platform'] = 'Illumina'
        elif any('ont' in f.lower() or 'nanopore' in f.lower() for f in file_matches):
            qc_data[sample_name]['platform'] = 'ONT'
        elif any('pacbio' in f.lower() for f in file_matches):
            qc_data[sample_name]['platform'] = 'PacBio'

        # Extract quality scores (Q30+ and Q40+)
        q_scores = re.findall(r'Q(\\d+): (\\d+)', content)
        q30_plus = sum(int(count) for q, count in q_scores if int(q) >= 30)
        q40_plus = sum(int(count) for q, count in q_scores if int(q) >= 40)

        qc_data[sample_name]['q30_plus'] = q30_plus
        qc_data[sample_name]['q40_plus'] = q40_plus

# Write TSV summary table
with open('qc_summary_table.tsv', 'w') as f:
    # Header
    f.write('Sample\\tPlatform\\tTotal_Reads\\tTotal_Bases\\tAvg_Read_Length\\tQ30_Plus_Bases\\tQ40_Plus_Bases\\tFiles\\n')

    # Data rows
    for sample_name in sorted(qc_data.keys()):
        data = qc_data[sample_name]
        files_str = ';'.join(data['files'])
        f.write(f"{data['sample']}\\t{data['platform']}\\t{data['total_reads']}\\t"
               f"{data['total_bases']}\\t{data['avg_length']:.2f}\\t{data['q30_plus']}\\t"
               f"{data['q40_plus']}\\t{files_str}\\n")

print(f"Processed {len(qc_data)} QC files")
print("Generated qc_summary_table.tsv")
EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    echo -e "Sample\\tPlatform\\tTotal_Reads\\tTotal_Bases\\tAvg_Read_Length\\tQ30_Plus_Bases\\tQ40_Plus_Bases\\tFiles" > qc_summary_table.tsv
    echo -e "test_sample\\tIllumina\\t1000\\t150000\\t150.0\\t800\\t600\\ttest_R1.fastq.gz;test_R2.fastq.gz" >> qc_summary_table.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
