# Pipeline Output Structure

### The pipeline outputs are organized into the following directory structure:
```
results/
├── genomes/
│   ├── sample1.fasta
│   ├── sample2.fasta
│   └── sample3.fasta
├── ont_reads/
│   ├── sample1_ont_0001.fastq.gz
│   ├── sample1_ont.maf.gz
│   ├── sample2_ont_0001.fastq.gz
│   └── sample2_ont.maf.gz
├── pacbio_reads/
│   ├── sample1_pacbio_0001.fastq.gz
│   ├── sample1_pacbio.maf.gz
│   ├── sample2_pacbio_0001.fastq.gz
│   └── sample2_pacbio.maf.gz
├── illumina_reads/
│   ├── sample1_illumina_R1.fastq.gz
│   ├── sample1_illumina_R2.fastq.gz
│   ├── sample2_illumina_R1.fastq.gz
│   └── sample2_illumina_R2.fastq.gz
├── multiqc/
│   ├── multiqc_report.html
│   ├── multiqc_data/
│   └── multiqc_plots/
└── pipeline_info/
    ├── execution_report.html
    ├── execution_timeline.html
    ├── execution_trace.txt
    └── pipeline_dag.svg
```
Output Files Description
Reference Genomes (`genomes/`)
File	Description
*.fasta	Reference genome sequences downloaded from NCBI or copied from local files

Example:

| sample1.fasta | Human reference genome (GRCh38) |
|---|---|
| ecoli.fasta | E. coli reference genome |

ONT Reads (`ont_reads/`)

Generated when ont_reads > 0 in the input samplesheet.
File Pattern	Description
*_ont_*.fastq.gz	Simulated ONT long reads in FASTQ format (compressed)
*_ont.maf.gz	Multiple Alignment Format file showing read alignments to reference (optional)

File Details:

    Format: FASTQ (gzipped)
    Read Length: Configurable (default: mean=10kb, sd=5kb)
    Accuracy: Configurable (default: 95%)
    Simulator: PBSIM3 or NanoSim

Example Files:
sample1_ont_0001.fastq.gz    # ONT reads for sample1
sample1_ont.maf.gz           # Alignment information
PacBio Reads (pacbio_reads/)

Generated when pacbio_reads > 0 in the input samplesheet.
File Pattern	Description
*_pacbio_*.fastq.gz	Simulated PacBio long reads in FASTQ format (compressed)
*_pacbio.maf.gz	Multiple Alignment Format file showing read alignments to reference (optional)

File Details:

    Format: FASTQ (gzipped)
    Read Length: Configurable (default: mean=8kb, sd=3kb)
    Accuracy: Configurable (default: 87%)
    Simulator: PBSIM3

Example Files:
sample1_pacbio_0001.fastq.gz    # PacBio reads for sample1
sample1_pacbio.maf.gz           # Alignment information
Illumina Reads (illumina_reads/)

Generated when illumina_reads > 0 in the input samplesheet.
File Pattern	Description
*_illumina_R1.fastq.gz	Forward reads (R1) in FASTQ format (compressed)
*_illumina_R2.fastq.gz	Reverse reads (R2) in FASTQ format (compressed)
*.aln	Alignment information (optional, if not using --noALN)

File Details:

    Format: Paired-end FASTQ (gzipped)
    Read Length: Configurable (default: 150bp)
    Insert Size: Configurable (default: mean=400bp, sd=50bp)
    Platform: Configurable (default: HiSeq 2500)
    Simulator: ART

Example Files:
sample1_illumina_R1.fastq.gz    # Forward reads
sample1_illumina_R2.fastq.gz    # Reverse reads
Quality Control Reports (multiqc/)
File	Description
multiqc_report.html	Comprehensive HTML report summarizing all pipeline outputs
multiqc_data/	Directory containing raw data used for the report
multiqc_plots/	Directory containing plot files (if generated)
Pipeline Information (pipeline_info/)
File	Description
execution_report.html	Nextflow execution report with resource usage and timing
execution_timeline.html	Timeline visualization of pipeline execution
execution_trace.txt	Detailed trace of all executed processes
pipeline_dag.svg	Directed Acyclic Graph showing pipeline workflow
File Formats
FASTQ Format

All simulated reads are provided in standard FASTQ format:
@read_identifier
SEQUENCE
+
QUALITY_SCORES

Example ONT read:
@S1_000001/1_10245_0_10245_0_0_0_0.000000_0
ATCGATCGATCGATCG...
+
################...

Example Illumina read pair:
@sample1_illumina_1
ATCGATCGATCGATCG...
+
IIIIIIIIIIIIIIII...

@sample1_illumina_2
CGATCGATCGATCGAT...
+
IIIIIIIIIIIIIIII...
MAF Format (Multiple Alignment Format)

MAF files contain alignment information between simulated reads and the reference genome:
a score=0
s ref 0 100 + 1000 ATCGATCG...
s read 0 100 + 100 ATCGATCG...
Read Statistics
Expected Read Characteristics
Platform	Read Length	Error Rate	Error Types
ONT	1-50kb (configurable)	~5% (configurable)	Insertions, deletions, substitutions
PacBio	1-30kb (configurable)	~13% (configurable)	Insertions, deletions, substitutions
Illumina	50-300bp (configurable)	~0.1%	Substitutions primarily
Quality Score Encoding

    ONT/PacBio: Phred+33 encoding
    Illumina: Phred+33 encoding (standard)

Data Usage Examples
Loading Reads in Python
import gzip
from Bio import SeqIO

# Load ONT reads
with gzip.open('ont_reads/sample1_ont_0001.fastq.gz', 'rt') as f:
    ont_reads = list(SeqIO.parse(f, 'fastq'))

# Load Illumina paired reads
with gzip.open('illumina_reads/sample1_illumina_R1.fastq.gz', 'rt') as f:
    r1_reads = list(SeqIO.parse(f, 'fastq'))

with gzip.open('illumina_reads/sample1_illumina_R2.fastq.gz', 'rt') as f:
    r2_reads = list(SeqIO.parse(f, 'fastq'))
Loading Reads in R
library(ShortRead)

# Load Illumina reads
r1 <- readFastq("illumina_reads/sample1_illumina_R1.fastq.gz")
r2 <- readFastq("illumina_reads/sample1_illumina_R2.fastq.gz")

# Load ONT reads
ont_reads <- readFastq("ont_reads/sample1_ont_0001.fastq.gz")
File Size Estimates

Approximate file sizes for different numbers of reads:
Platform	Reads	Uncompressed Size	Compressed Size
ONT	10,000	~400 MB	~100 MB
PacBio	10,000	~320 MB	~80 MB
Illumina	1,000,000 pairs	~600 MB	~150 MB

Note: Actual sizes depend on read length, quality scores, and genome complexity.
Troubleshooting
Common Issues

    Empty output files: Check that the input parameters specify non-zero read counts
    Missing MAF files: MAF files are optional and may not be generated depending on simulator settings
    Unexpected file sizes: Verify read count parameters and expected read lengths

Validation

To validate your output:
# Check read counts
zcat ont_reads/*.fastq.gz | grep -c "^@"
zcat illumina_reads/*_R1.fastq.gz | grep -c "^@"

# Check read length distribution
zcat ont_reads/*.fastq.gz | awk 'NR%4==2{print length($0)}' | sort -n | uniq -c
Next Steps

The simulated reads can be used for:

    Pipeline testing: Validate bioinformatics workflows
    Method benchmarking: Compare analysis tools and parameters
    Algorithm development: Test new computational methods
    Training datasets: Educational purposes and tutorials

For downstream analysis, consider using established pipelines such as:

    nf-core/nanoseq for ONT/PacBio analysis
    nf-core/rnaseq for RNA-seq analysis
    nf-core/sarek for variant calling

For questions about output interpretation or pipeline usage, please refer to the main documentation or open an issue on the project repository.
