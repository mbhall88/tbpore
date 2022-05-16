# TBpore

*Mycobacterium tuberculosis* genomic analysis from Nanopore sequencing data

[![Python CI](https://github.com/mbhall88/tbpore/actions/workflows/ci.yaml/badge.svg)](https://github.com/mbhall88/tbpore/actions/workflows/ci.yaml)
[![codecov](https://codecov.io/gh/mbhall88/tbpore/branch/main/graph/badge.svg)](https://codecov.io/gh/mbhall88/tbpore)
[![PyPI](https://img.shields.io/pypi/v/tbpore)](https://pypi.org/project/tbpore/)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/tbpore)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

**⚠️WORK IN PROGRESS⚠️**

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
scripts/run_sample_example.sh
```

# Performance

Benchmarked on 91 TB Madagascar ONT samples with 1 thread:
* Runtime: `2255`s avg, `3805`s max (s = seconds);
* RAM: `12.4`GB avg, `13.1`GB max (GB = Gigabytes);

# Usage

```
Usage: tbpore [OPTIONS] [INPUTS]...

  Mycobacterium tuberculosis genomic analysis from Nanopore sequencing data

  INPUTS: Fastq file(s) and/or a directory containing fastq files. All files
  will be joined into a single fastq file, so ensure they're all part of the
  same sample/isolate.

Options:
  -h, --help                      Show this message and exit.
  -V, --version                   Show the version and exit.
  -o, --outdir DIRECTORY          Directory to place output files  [default:
                                  tbpore_out]
  -v, --verbose                   Turns on debug-level logger. Option is
                                  mutually exclusive with quiet.
  -q, --quiet                     Turns off all logging except errors. Option
                                  is mutually exclusive with verbose.
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
```