process CREATE_MANIFEST {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(ont_fastq), path(pacbio_fastq), path(illumina_r1), path(illumina_r2)

    output:
    tuple val(meta), path("${meta.id}_manifest.csv"), emit: manifest
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Handle optional file inputs - check if files exist and are not empty placeholders
    def ont_path = ont_fastq && ont_fastq.name != 'NO_FILE' && ont_fastq.size() > 0 ? ont_fastq : ''
    def pacbio_path = pacbio_fastq && pacbio_fastq.name != 'NO_FILE' && pacbio_fastq.size() > 0 ? pacbio_fastq : ''
    def illumina_r1_path = illumina_r1 && illumina_r1.name != 'NO_FILE' && illumina_r1.size() > 0 ? illumina_r1 : ''
    def illumina_r2_path = illumina_r2 && illumina_r2.name != 'NO_FILE' && illumina_r2.size() > 0 ? illumina_r2 : ''

    """
    #!/usr/bin/env python3

    import csv
    import os

    # Define the sample name
    sample_name = "${prefix}"

    # Define file paths (empty string if not provided)
    ont_fastq = "${ont_path}"
    pacbio_fastq = "${pacbio_path}"
    illumina_r1 = "${illumina_r1_path}"
    illumina_r2 = "${illumina_r2_path}"

    # Create the header dynamically based on which files are present
    header = ["sample"]
    data_row = [sample_name]

    # Add columns only if the corresponding files exist
    if ont_fastq:
        header.append("ont_fastq")
        data_row.append(ont_fastq)

    if pacbio_fastq:
        header.append("pacbio_fastq")
        data_row.append(pacbio_fastq)

    if illumina_r1:
        header.append("illumina_fastq_R1")
        data_row.append(illumina_r1)

    if illumina_r2:
        header.append("illumina_fastq_R2")
        data_row.append(illumina_r2)

    # Write the manifest CSV
    with open("${prefix}_manifest.csv", "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)
        writer.writerow(data_row)

    # Create versions file
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write('    python: "3.9"\\n')
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_manifest.csv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.9"
    END_VERSIONS
    """
}
