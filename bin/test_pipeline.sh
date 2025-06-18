#!/bin/bash

# SeqAudit Pipeline Test Script
# This script runs the pipeline test profiles to validate the installation

set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    local color=$1
    local message=$2
    echo -e "${color}${message}${NC}"
}

# Function to run a test
run_test() {
    local profile=$1
    local description=$2
    local outdir="test_results_${profile}"
    
    print_status $YELLOW "Running ${description}..."
    
    if nextflow run . -profile ${profile},docker --outdir ${outdir} -resume; then
        print_status $GREEN "✓ ${description} completed successfully"
        return 0
    else
        print_status $RED "✗ ${description} failed"
        return 1
    fi
}

# Main execution
main() {
    print_status $YELLOW "SeqAudit Pipeline Test Suite"
    print_status $YELLOW "============================"
    echo
    
    # Check if Nextflow is available
    if ! command -v nextflow &> /dev/null; then
        print_status $RED "Error: Nextflow is not installed or not in PATH"
        exit 1
    fi
    
    # Check if Docker is available
    if ! command -v docker &> /dev/null; then
        print_status $RED "Error: Docker is not installed or not in PATH"
        exit 1
    fi
    
    # Validate test samplesheets
    print_status $YELLOW "Validating test samplesheets..."
    for samplesheet in assets/samplesheet_test*.csv; do
        if [[ -f "$samplesheet" ]]; then
            if python3 bin/validate_samplesheet.py "$samplesheet"; then
                print_status $GREEN "✓ $(basename $samplesheet) is valid"
            else
                print_status $RED "✗ $(basename $samplesheet) validation failed"
                exit 1
            fi
        fi
    done
    echo
    
    # Run tests
    local failed_tests=0
    
    # CI Test (fastest)
    if ! run_test "test_ci" "CI Test (minimal resources)"; then
        ((failed_tests++))
    fi
    echo
    
    # Basic Test
    if ! run_test "test" "Basic Test"; then
        ((failed_tests++))
    fi
    echo
    
    # Ask user if they want to run full test
    read -p "Run full test? This may take several hours (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        if ! run_test "test_full" "Full Test"; then
            ((failed_tests++))
        fi
        echo
    fi
    
    # Summary
    print_status $YELLOW "Test Summary"
    print_status $YELLOW "============"
    
    if [[ $failed_tests -eq 0 ]]; then
        print_status $GREEN "All tests passed! ✓"
        print_status $GREEN "SeqAudit pipeline is ready to use."
    else
        print_status $RED "${failed_tests} test(s) failed ✗"
        print_status $RED "Please check the error messages above."
        exit 1
    fi
}

# Show help
show_help() {
    echo "SeqAudit Pipeline Test Script"
    echo
    echo "Usage: $0 [OPTIONS]"
    echo
    echo "Options:"
    echo "  -h, --help     Show this help message"
    echo "  --ci-only      Run only CI test"
    echo "  --basic-only   Run only basic test"
    echo "  --full-only    Run only full test"
    echo
    echo "Examples:"
    echo "  $0                # Run interactive test suite"
    echo "  $0 --ci-only      # Run only CI test"
    echo "  $0 --basic-only   # Run only basic test"
}

# Parse command line arguments
case "${1:-}" in
    -h|--help)
        show_help
        exit 0
        ;;
    --ci-only)
        run_test "test_ci" "CI Test (minimal resources)"
        ;;
    --basic-only)
        run_test "test" "Basic Test"
        ;;
    --full-only)
        run_test "test_full" "Full Test"
        ;;
    "")
        main
        ;;
    *)
        print_status $RED "Unknown option: $1"
        show_help
        exit 1
        ;;
esac