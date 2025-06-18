#!/usr/bin/env python3

"""
Validate SeqAudit samplesheet format
"""

import argparse
import csv
import sys
from pathlib import Path


def validate_samplesheet(samplesheet_path):
    """Validate the samplesheet format and content"""
    
    required_columns = ['sample_id', 'genome_source', 'genome_id', 'ont_reads', 'pacbio_reads', 'illumina_reads']
    valid_genome_sources = ['refseq', 'genbank', 'local']
    
    errors = []
    warnings = []
    
    try:
        with open(samplesheet_path, 'r') as f:
            reader = csv.DictReader(f)
            
            # Check headers
            if not all(col in reader.fieldnames for col in required_columns):
                missing = [col for col in required_columns if col not in reader.fieldnames]
                errors.append(f"Missing required columns: {', '.join(missing)}")
                return errors, warnings
            
            # Check each row
            for i, row in enumerate(reader, start=2):  # Start at 2 because header is row 1
                
                # Check sample_id
                if not row['sample_id'].strip():
                    errors.append(f"Row {i}: sample_id cannot be empty")
                
                # Check genome_source
                if row['genome_source'] not in valid_genome_sources:
                    errors.append(f"Row {i}: genome_source must be one of {valid_genome_sources}")
                
                # Check genome_id
                if not row['genome_id'].strip():
                    errors.append(f"Row {i}: genome_id cannot be empty")
                
                # Check read counts are integers
                for read_type in ['ont_reads', 'pacbio_reads', 'illumina_reads']:
                    try:
                        count = int(row[read_type])
                        if count < 0:
                            errors.append(f"Row {i}: {read_type} must be >= 0")
                    except ValueError:
                        errors.append(f"Row {i}: {read_type} must be an integer")
                
                # Check if local file exists
                if row['genome_source'] == 'local':
                    if not Path(row['genome_id']).exists():
                        warnings.append(f"Row {i}: Local file {row['genome_id']} does not exist")
                
                # Check if at least one read type is requested
                try:
                    ont = int(row['ont_reads'])
                    pacbio = int(row['pacbio_reads'])
                    illumina = int(row['illumina_reads'])
                    if ont + pacbio + illumina == 0:
                        warnings.append(f"Row {i}: No reads requested for sample {row['sample_id']}")
                except ValueError:
                    pass  # Already caught above
    
    except FileNotFoundError:
        errors.append(f"Samplesheet file not found: {samplesheet_path}")
    except Exception as e:
        errors.append(f"Error reading samplesheet: {str(e)}")
    
    return errors, warnings


def main():
    parser = argparse.ArgumentParser(description='Validate SeqAudit samplesheet')
    parser.add_argument('samplesheet', help='Path to samplesheet CSV file')
    parser.add_argument('--strict', action='store_true', help='Treat warnings as errors')
    
    args = parser.parse_args()
    
    errors, warnings = validate_samplesheet(args.samplesheet)
    
    if errors:
        print("ERRORS:", file=sys.stderr)
        for error in errors:
            print(f"  - {error}", file=sys.stderr)
    
    if warnings:
        print("WARNINGS:", file=sys.stderr)
        for warning in warnings:
            print(f"  - {warning}", file=sys.stderr)
    
    if errors or (args.strict and warnings):
        sys.exit(1)
    else:
        print("Samplesheet validation passed!")
        sys.exit(0)


if __name__ == '__main__':
    main()