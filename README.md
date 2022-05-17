# TBpore

*Mycobacterium tuberculosis* genomic analysis from Nanopore sequencing data

[![Python CI](https://github.com/mbhall88/tbpore/actions/workflows/ci.yaml/badge.svg)](https://github.com/mbhall88/tbpore/actions/workflows/ci.yaml)
[![codecov](https://codecov.io/gh/mbhall88/tbpore/branch/main/graph/badge.svg)](https://codecov.io/gh/mbhall88/tbpore)
[![PyPI](https://img.shields.io/pypi/v/tbpore)](https://pypi.org/project/tbpore/)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/tbpore)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

**⚠️WORK IN PROGRESS⚠️**

[TOC]: #

# Table of Contents
- [Synopsis](#synopsis)
- [Installation](#installation)
- [Configuring the decontamination database index](#configuring-the-decontamination-database-index)
- [Performance](#performance)
- [Usage](#usage)

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

<!---
## conda

`conda install tbpore`
-->

## local

### Dependencies
* `git`
* `conda`

### Walkthrough

The following steps will create a [conda](https://docs.conda.io/en/latest/) environment named `tbpore` which will
contain all dependencies and `tbpore` itself.

```
git clone https://github.com/mbhall88/tbpore  # get tbpore source code
cd tbpore
conda env create -f environment.yaml && conda activate tbpore  # install dependencies
just install  # install tbpore
just check  # check installation is fine
tbpore -h  # print usage
```

Whenever you want to rerun `tbpore` and you are not already in the `conda` `tbpore` environment, you can activate it by
running `conda activate tbpore` and then `tbpore` will be available.

# Configuring the decontamination database index

When you run your first `tbpore process`, you will get this error:
```
ERROR    | Decontamination DB index tbpore/data/decontamination_db/tbpore.remove_contam.fa.gz.map-ont.mmi does not
exist, please follow the instructions at https://github.com/mbhall88/tbpore#configuring-the-decontamination-database-index
to download and configure it before running tbpore
```
This means you need to download the [minimap2](https://github.com/lh3/minimap2) decontamination database index before
proceeding. We did not include this index in this repo as it is too heavy. Right now this index is private and can be
downloaded by logging with your credentials into the EBI private FTP server `ftp-private.ebi.ac.uk` and downloading the
index at `tbpore/0.1.0/decontamination_db/tbpore.remove_contam.fa.gz.map-ont.mmi.gz`. Once the download is complete,
you can:

1. Ensure that the compressed index was transferred correctly by checking its `md5sum`:
```
$ md5sum tbpore.remove_contam.fa.gz.map-ont.mmi.gz
82d050e0f1cba052f0c94f16fcb32f7b  tbpore.remove_contam.fa.gz.map-ont.mmi.gz
```

2. Decompress the index:
```
gunzip tbpore.remove_contam.fa.gz.map-ont.mmi.gz
```

3. Check the md5sum of the decompressed index:
```
$ md5sum tbpore.remove_contam.fa.gz.map-ont.mmi
810c5c09eaf9421128e4e52cdf2fa32a  tbpore.remove_contam.fa.gz.map-ont.mmi
```

4. Move the decompressed index to `<tbpore_dir>/data/decontamination_db/tbpore.remove_contam.fa.gz.map-ont.mmi`

Once these four steps above are done, you should be able to run `tbpore` on an example isolate by going into the
`tbpore` dir and running:
```
just test-run
```

# Performance

## `tbpore process`

Benchmarked on 151 TB ONT samples with 1 thread:
* Runtime: `2103`s avg, `4048`s max (s = seconds);
* RAM: `12.4`GB avg, `13.1`GB max (GB = Gigabytes);

## `tbpore cluster`

Clustering 151 TB ONT samples:
* Runtime: `286`s;
* RAM: `<1`GB;

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