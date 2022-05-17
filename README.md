# TBpore

*Mycobacterium tuberculosis* genomic analysis from Nanopore sequencing data

[![Python CI](https://github.com/mbhall88/tbpore/actions/workflows/ci.yaml/badge.svg)](https://github.com/mbhall88/tbpore/actions/workflows/ci.yaml)
[![codecov](https://codecov.io/gh/mbhall88/tbpore/branch/main/graph/badge.svg)](https://codecov.io/gh/mbhall88/tbpore)
[![PyPI](https://img.shields.io/pypi/v/tbpore)](https://pypi.org/project/tbpore/)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/tbpore)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

**⚠️WORK IN PROGRESS⚠️**

# Synopsis

`tbpore` is a tool with two main goals.
First is to process Nanopore Mycobacterium tuberculosis sequencing data to describe variants with respect to the
canonical TB strain H37Rv and predict antibiotic resistance (command `tbpore process`).
Variant description is done by decontaminating reads, calling variants with
[bcftools](https://github.com/samtools/bcftools) and filtering variants.
Antibiotic resistance is predicted with [mykrobe](https://github.com/Mykrobe-tools/mykrobe).
Second, `tbpore` can be used to cluster TB samples based on their genotyping and a given distance threshold (command
`tbpore cluster`).

# Installation

## conda

`conda install tbpore`

## local

### Dependencies
* `git`
* `conda`

### Walkthrough

```
git clone https://github.com/mbhall88/tbpore
cd tbpore
conda env create -f environment.yaml && conda activate tbpore  # install dependencies
just install  # install tbpore
just check  # checks installation is fine
# if you want to run tbpore on an example isolate (will require you to download the decontamination DB index)
just test-run
```

# Performance

Benchmarked on 91 TB Madagascar ONT samples with 1 thread:
* Runtime: `2255`s avg, `3805`s max (s = seconds);
* RAM: `12.4`GB avg, `13.1`GB max (GB = Gigabytes);

# Usage

## General usage

```
Usage: tbpore [OPTIONS] COMMAND [ARGS]...

Options:
  -h, --help     Show this message and exit.
  -V, --version  Show the version and exit.
  -v, --verbose  Turns on debug-level logger. Option is mutually exclusive
                 with quiet.
  -q, --quiet    Turns off all logging except errors. Option is mutually
                 exclusive with verbose.

Commands:
  cluster  Cluster consensus sequences
  process  Single-sample TB genomic analysis from Nanopore sequencing data
```

## process subcommand

```
Usage: tbpore process [OPTIONS] [INPUTS]...

  Single-sample TB genomic analysis from Nanopore sequencing data

  INPUTS: Fastq file(s) and/or a directory containing fastq files. All files
  will be joined into a single fastq file, so ensure they're all part of the
  same sample/isolate.

Options:
  -o, --outdir DIRECTORY          Directory to place output files  [default:
                                  tbpore_out]
  -r, --recursive                 Recursively search INPUTS for fastq files
  --tmp DIRECTORY                 Specify where to write all (tbpore)
                                  temporary files. [default: <outdir>/.tbpore]
  -S, --name TEXT                 Name of the sample. By default, will use the
                                  first INPUT file with any extensions
                                  stripped
  -t, --threads INTEGER           Number of threads to use in multithreaded
                                  tools  [default: 1]
  -A, --report_all_mykrobe_calls  Report all mykrobe calls (turn on flag -A,
                                  --report_all_calls when calling mykrobe)
  -d, --cleanup / -D, --no-cleanup
                                  Remove all temporary files on *successful*
                                  completion  [default: no-cleanup]
  --help                          Show this message and exit.
```

## cluster subcommand

```
Usage: tbpore cluster [OPTIONS] [INPUTS]...

  Cluster consensus sequences

  Preferably input consensus sequences previously generated with tbpore
  process.

  INPUTS: Two or more consensus fasta sequences. Use glob patterns to input
  several easily (e.g. output/sample_*/*.consensus.fa).

Options:
  -T, --threshold INTEGER         Clustering threshold  [default: 6]
  -o, --outdir DIRECTORY          Directory to place output files  [default:
                                  cluster_out]
  --tmp DIRECTORY                 Specify where to write all (tbpore)
                                  temporary files. [default: <outdir>/.tbpore]
  -t, --threads INTEGER           Number of threads to use in multithreaded
                                  tools  [default: 1]
  -d, --cleanup / -D, --no-cleanup
                                  Remove all temporary files on *successful*
                                  completion  [default: no-cleanup]
  --help                          Show this message and exit.
```