# nf-core/seqaudit

[![GitHub Actions CI Status](https://github.com/nf-core/seqaudit/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-core/seqaudit/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/seqaudit/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/seqaudit/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/seqaudit)

## Introduction

**SeqAudit** is a bioinformatics pipeline that simulating ONT, PacBio, and Illumina sequencing data'

Create a set of fastq files from an input reference genome fasta. You can provide the number of sample sets that you wish to create for each technology type (Illumina fastq will be created as paired end reads) and provide the reference genome as a path, RefSeq, or Genbank ID,

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:
##### Example samplesheet.csv
```
sample_id,genome_source,genome_id,ont_reads,pacbio_reads,illumina_reads
sample1,refseq,GCF_000001405.39,10000,5000,1000000
sample2,genbank,GCA_000005825.2,5000,8000,500000
sample3,local,/path/to/genome.fasta,0,10000,2000000
ecoli,refseq,GCF_000005825.2,15000,12000,800000
```
Each row represents the creation of X number of a fastq files (single-end) or pairs of fastq files (paired end) for a specific sequencing technology

Now, you can run the basic pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
# Simulate ONT, PacBio, and Illumina reads
nextflow run main.nf \
    --input samplesheet.csv \
    --outdir results \
    -profile docker
```

You can also provide more parameters for use with specific sequencing technologies
> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Credits

SeqAudit was originally written by Holly Halstead.



<!-- TODO nf-core: If applicable, make list of people who have also contributed
We thank the following people for their extensive assistance in the development of this pipeline:
-->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf-core/seqaudit for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
