# TBpore

*Mycobacterium tuberculosis* genomic analysis from Nanopore sequencing data

[![Python CI](https://github.com/mbhall88/tbpore/actions/workflows/ci.yaml/badge.svg)](https://github.com/mbhall88/tbpore/actions/workflows/ci.yaml)
[![codecov](https://codecov.io/gh/mbhall88/tbpore/branch/main/graph/badge.svg)](https://codecov.io/gh/mbhall88/tbpore)
[![PyPI](https://img.shields.io/pypi/v/tbpore)](https://pypi.org/project/tbpore/)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/tbpore)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

[TOC]: #

# Table of Contents
- [Synopsis](#synopsis)
- [Installation](#installation)
- [Walkthrough](#walkthrough)
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

## Installation

### conda

[![Conda (channel only)](https://img.shields.io/conda/vn/bioconda/tbpore)](https://anaconda.org/bioconda/tbpore)
[![bioconda version](https://anaconda.org/bioconda/tbpore/badges/platforms.svg)](https://anaconda.org/bioconda/tbpore)
![Conda](https://img.shields.io/conda/dn/bioconda/tbpore)

Prerequisite: [`conda`][conda] (and bioconda channel [correctly set up][channels])

```shell
$ conda install tbpore
```

### pip

![PyPI](https://img.shields.io/pypi/v/tbpore)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/tbpore)

The python components of `tbpore` are availble to install through [PyPI].

```shell
pip install tbpore
```

**However**, you will need to install the following dependencies, which cannot be installed through PyPI.

#### Dependencies
* [`rasusa`](https://github.com/mbhall88/rasusa)
* [`psdm`](https://github.com/mbhall88/psdm) version 0.1
* [`samtools`](https://github.com/samtools/samtools) version 1.13
* [`bcftools`](https://github.com/samtools/bcftools) version 1.13
* [`mykrobe`](https://github.com/Mykrobe-tools/mykrobe) version ≥ 0.12
* [`minimap2`](https://github.com/lh3/minimap2) version 2.22
* [`seqkit`](https://bioinf.shenwei.me/seqkit/) version 2.0

We make no guarentees about the performance of `tbpore` with versions other than those specified above. In particular, the `bcftools` version is very important. The latest versions of the other dependencies can likely be used.

### Container

Docker images are provided through biocontainers.

#### `singularity`

Prerequisite: [`singularity`][singularity]

```shell
$ URI="docker://quay.io/biocontainers/tbpore:<tag>"
$ singularity exec "$URI" tbpore --help
```

see [here][tags] for valid values for `<tag>`.

#### `docker`

[![Docker Repository on Quay](https://quay.io/repository/biocontainers/tbpore/status "Docker Repository on Quay")](https://quay.io/repository/biocontainers/tbpore)

Prerequisite: [Docker]

```shell
$ docker pull quay.io/biocontainers/tbpore:<tag>
$ docker run quay.io/biocontainers/tbpore:<tag> tbpore --help
```

see [here][tags] for valid values for `<tag>`.

### Configuring the decontamination database index

When you run your first `tbpore process`, you will get this error:
```
ERROR    | Decontamination DB index tbpore/data/decontamination_db/tbpore.remove_contam.fa.gz.map-ont.mmi does not
exist, please follow the instructions at https://github.com/mbhall88/tbpore#configuring-the-decontamination-database-index
to download and configure it before running tbpore
```
This means you need to download the [minimap2](https://github.com/lh3/minimap2) decontamination database index before
proceeding. You can [download this index here](https://figshare.com/ndownloader/files/36708444) or by running:
```shell
wget https://figshare.com/ndownloader/files/36708444 -O tbpore.remove_contam.fa.gz.map-ont.mmi.gz
```

Once the download is complete, you can:

1. Ensure that the compressed index was transferred correctly by checking its `md5sum`:
```shell
md5sum tbpore.remove_contam.fa.gz.map-ont.mmi.gz
82d050e0f1cba052f0c94f16fcb32f7b  tbpore.remove_contam.fa.gz.map-ont.mmi.gz
```

2. Decompress the index:
```shell
gunzip tbpore.remove_contam.fa.gz.map-ont.mmi.gz
```

3. Check the md5sum of the decompressed index:
```shell
md5sum tbpore.remove_contam.fa.gz.map-ont.mmi
810c5c09eaf9421128e4e52cdf2fa32a  tbpore.remove_contam.fa.gz.map-ont.mmi
```

4. Move the decompressed index to `<tbpore_dir>/data/decontamination_db/tbpore.remove_contam.fa.gz.map-ont.mmi`
    * Note: you can also keep this index at a different path and specify it to `tbpore` using the `--db` option;

Once these four steps above are done, you should be able to run `tbpore` on an example isolate by going into the
`tbpore` dir and running:
```shell
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
  -h, --help                      Show this message and exit.
  -r, --recursive                 Recursively search INPUTS for fastq files
  -S, --name TEXT                 Name of the sample. By default, will use the
                                  first INPUT file with fastq extensions
                                  removed
  -A, --report_all_mykrobe_calls  Report all mykrobe calls (turn on flag -A,
                                  --report_all_calls when calling mykrobe)
  --db PATH                       Path to the decontaminaton database
                                  [default: /Users/michaelhall/Projects/tbpore
                                  /data/decontamination_db/tbpore.remove_conta
                                  m.fa.gz.map-ont.mmi]
  -m, --metadata PATH             Path to the decontaminaton database metadata
                                  file  [default: /Users/michaelhall/Projects/
                                  tbpore/data/decontamination_db/remove_contam
                                  .tsv.gz]
  -o, --outdir DIRECTORY          Directory to place output files  [default:
                                  .]
  --tmp DIRECTORY                 Specify where to write all (tbpore)
                                  temporary files. [default: <outdir>/.tbpore]
  -t, --threads INTEGER           Number of threads to use in multithreaded
                                  tools  [default: 1]
  -d, --cleanup / -D, --no-cleanup
                                  Remove all temporary files on *successful*
                                  completion  [default: no-cleanup]
  --cache DIRECTORY               Path to use for the cache  [default:
                                  /Users/michaelhall/.cache]
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
  -h, --help                      Show this message and exit.
  -T, --threshold INTEGER         Clustering threshold  [default: 6]
  -o, --outdir DIRECTORY          Directory to place output files  [default:
                                  .]
  --tmp DIRECTORY                 Specify where to write all (tbpore)
                                  temporary files. [default: <outdir>/.tbpore]
  -t, --threads INTEGER           Number of threads to use in multithreaded
                                  tools  [default: 1]
  -d, --cleanup / -D, --no-cleanup
                                  Remove all temporary files on *successful*
                                  completion  [default: no-cleanup]
  --cache DIRECTORY               Path to use for the cache  [default:
                                  /Users/michaelhall/.cache]
```

[channels]: https://bioconda.github.io/#usage
[conda]: https://docs.conda.io/projects/conda/en/latest/user-guide/install/
[PyPI]: https://pypi.org/project/tbpore/
[singularity]: https://sylabs.io/guides/3.6/user-guide/quick_start.html#quick-installation-steps
[tags]: https://quay.io/repository/biocontainers/tbpore?tab=tags
[Docker]: https://docs.docker.com/v17.12/install/