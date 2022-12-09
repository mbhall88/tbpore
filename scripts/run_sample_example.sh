#!/usr/bin/env sh
set -eu

reads_url="ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR903/ERR9030361/mada_1-7.subsampled.fastq.gz"

mkdir -p sample_example

reads="sample_example/mada_1-7.subsampled.fastq.gz"
if [ ! -f "$reads" ]; then
    wget "$reads_url" -O "$reads"
fi

poetry run tbpore process -o sample_example/tbpore_out --cleanup "$reads"
