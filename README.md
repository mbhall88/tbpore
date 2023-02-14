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
- [Citation](#citation)
- [Installation](#installation)
- [Configuring the decontamination database index](#configuring-the-decontamination-database-index)
- [Performance](#performance)
- [Usage](#usage)

# Synopsis

`tbpore` is a tool with two main goals.
First is to process Nanopore Mycobacterium tuberculosis sequencing data to describe
variants with respect to the
canonical TB strain H37Rv and predict antibiotic resistance (command `tbpore process`).
Variant description is done by decontaminating reads, calling variants with
[bcftools](https://github.com/samtools/bcftools) and filtering variants.
Antibiotic resistance is predicted
with [mykrobe](https://github.com/Mykrobe-tools/mykrobe).
Second, `tbpore` can be used to cluster TB samples based on their genotyping and a given
distance threshold (command
`tbpore cluster`).

## Citation

TBpore is a slimmed-down version of
the [full pipeline](https://github.com/mbhall88/head_to_head_pipeline) used
in our paper 👇


> Hall, M. B. et al. Evaluation of Nanopore sequencing for Mycobacterium tuberculosis drug susceptibility testing and outbreak investigation: a genomic analysis. *The Lancet Microbe* 0, (2022) doi: [10.1016/S2666-5247(22)00301-9][doi].

[doi]: https://doi.org/10.1016/S2666-5247(22)00301-9

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

**However**, you will need to install the following dependencies, which cannot be
installed through PyPI.

#### Dependencies

* [`rasusa`](https://github.com/mbhall88/rasusa)
* [`psdm`](https://github.com/mbhall88/psdm) version 0.1.x
* [`samtools`](https://github.com/samtools/samtools) version 1.13
* [`bcftools`](https://github.com/samtools/bcftools) version 1.13
* [`mykrobe`](https://github.com/Mykrobe-tools/mykrobe) version 0.12.x
* [`minimap2`](https://github.com/lh3/minimap2) version 2.22
* [`seqkit`](https://bioinf.shenwei.me/seqkit/) version 2.x
* [`nanoq`](https://github.com/esteinig/nanoq) version 0.9.x

We make no guarentees about the performance of `tbpore` with versions other than those
specified above. In particular, the `bcftools` version is very important. The latest
versions of the other dependencies can likely be used.

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

## Configuring the decontamination database index

After installing TBpore, you will need to download the decontamination database index.

```
$ tbpore download
```

By default, this will download the index
to `${HOME}/.tbpore/decontamination_db/remove_contam.map-ont.mmi`, as this is the
default location `tbpore process` will search for.

If you prefer to download the index to another location, this can be done with

```
$ tbpore download -o other/location/db.mmi
```

Keep in mind, if you specify a non-default location, you will need to use the `--db`
option when running `tbpore process`.

## Performance

### `tbpore process`

Benchmarked on 151 TB ONT samples with 1 thread:

* Runtime: `2103`s avg, `4048`s max (s = seconds);
* RAM: `12.4`GB avg, `13.1`GB max (GB = Gigabytes);

### `tbpore cluster`

Clustering 151 TB ONT samples:

* Runtime: `286`s;
* RAM: `<1`GB;

## Usage

### General usage

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
  cluster   Cluster consensus sequences
  download  Download and validate the decontamination database
  process   Single-sample TB genomic analysis from Nanopore sequencing data
```

### process

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
                                  [default: ${HOME}/.tbpore/decontamination_db/
                                  remove_contam.map-ont.mmi]
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

### cluster

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

### download

```
Usage: tbpore download [OPTIONS]

  Download and validate the decontamination database

Options:
  -h, --help         Show this message and exit.
  -o, --output PATH  Download database to a specified filepath  [default: ${HOME}/
                     .tbpore/decontamination_db/remove_contam.map-ont.mmi]
  -f, --force        Force overwrite if the database already exists
```

[channels]: https://bioconda.github.io/#usage

[conda]: https://docs.conda.io/projects/conda/en/latest/user-guide/install/

[PyPI]: https://pypi.org/project/tbpore/

[singularity]: https://sylabs.io/guides/3.6/user-guide/quick_start.html#quick-installation-steps

[tags]: https://quay.io/repository/biocontainers/tbpore?tab=tags

[Docker]: https://docs.docker.com/v17.12/install/
