#!/bin/env sh
set -e
set -u

reads="ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR903/ERR9030361/mada_1-7.subsampled.fastq.gz"

mkdir -p sample_example
wget $reads -O sample_example/mada_1-7.subsampled.fastq.gz
tbpore -o sample_example/tbpore_out --cleanup sample_example/mada_1-7.subsampled.fastq.gz
