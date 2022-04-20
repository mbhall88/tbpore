import click
from loguru import logger

from tbpore import __version__, TMP_NAME, H37RV_genome
from tbpore.utils import concatenate_fastqs, find_fastq_files
from tbpore.external_tools import ExternalTool

import gzip
import shutil
import sys
from pathlib import Path
from typing import Tuple
import subprocess


log_fmt = (
    "[<green>{time:YYYY-MM-DD HH:mm:ss}</green>] <level>{level: <8}</level> | "
    "<level>{message}</level>"
)


class Mutex(click.Option):
    def __init__(self, *args, **kwargs):
        self.not_required_if: list = kwargs.pop("not_required_if")

        assert self.not_required_if, "'not_required_if' parameter required"
        kwargs["help"] = (
            kwargs.get("help", "")
            + "Option is mutually exclusive with "
            + ", ".join(self.not_required_if)
            + "."
        ).strip()
        super(Mutex, self).__init__(*args, **kwargs)

    def handle_parse_result(self, ctx, opts, args):
        current_opt: bool = self.name in opts
        for mutex_opt in self.not_required_if:
            if mutex_opt in opts:
                if current_opt:
                    raise click.UsageError(
                        "Illegal usage: '--"
                        + str(self.name)
                        + "' is mutually exclusive with '--"
                        + str(mutex_opt)
                        + "'"
                    )
                else:
                    self.prompt = None
        return super(Mutex, self).handle_parse_result(ctx, opts, args)


@click.command()
@click.help_option("--help", "-h")
@click.version_option(__version__, "--version", "-V")
@click.option(
    "-o",
    "--outdir",
    help="Directory to place output files",
    default=".",
    show_default=True,
    type=click.Path(file_okay=False, writable=True, path_type=Path),
)
@click.option(
    "-v",
    "--verbose",
    help="Turns on debug-level logger. ",
    is_flag=True,
    cls=Mutex,
    not_required_if=["quiet"],
)
@click.option(
    "-q",
    "--quiet",
    help="Turns off all logging except errors. ",
    is_flag=True,
    cls=Mutex,
    not_required_if=["verbose"],
)
@click.option(
    "-r", "--recursive", help="Recursively search INPUTS for fastq files", is_flag=True
)
@click.option(
    "--tmp",
    help=(
        f"Specify where to write all (tbpore) temporary files. [default: "
        f"<outdir>/{TMP_NAME}]"
    ),
    type=click.Path(file_okay=False, writable=True, path_type=Path),
)
@click.option(
    "-S",
    "--name",
    help=(
        "Name of the sample. By default, will use the first INPUT file with any "
        "extensions stripped"
    ),
)
@click.option(
    "-t",
    "--threads",
    help="Number of threads to use in multithreaded tools",
    type=int,
    default=1
)
@click.option(
    "-m",
    "--mem_mb",
    help="Memory to use in tools that accept a memory limit",
    type=int,
    default=8*1024  # 8 GB
)
@click.option(
    "--cleanup/--no-cleanup",
    "-d/-D",
    default=False,
    show_default=True,
    help="Remove all temporary files on *successful* completion",
)
@click.argument("inputs", type=click.Path(exists=True, path_type=Path), nargs=-1)
@click.pass_context
def main(
    ctx: click.Context,
    verbose: bool,
    quiet: bool,
    outdir: Path,
    inputs: Tuple[Path, ...],
    recursive: bool,
    tmp: Path,
    name: str,
    threads: int,
    mem_mb: int,
    cleanup: bool,
):
    """Mycobacterium tuberculosis genomic analysis from Nanopore sequencing data

    INPUTS: Fastq file(s) and/or a directory containing fastq files. All files will
    be joined into a single fastq file, so ensure thery're all part of the same
    sample/isolate.
    """
    log_lvl = "INFO"
    if verbose:
        log_lvl = "DEBUG"
    elif quiet:
        log_lvl = "ERROR"
    logger.remove()
    logger.add(sys.stderr, level=log_lvl, format=log_fmt)
    logger.info(f"Welcome to TBpore version {__version__}")

    outdir.mkdir(exist_ok=True, parents=True)
    if tmp is None:
        tmp = outdir / TMP_NAME
    tmp.mkdir(exist_ok=True, parents=True)

    if not inputs:
        logger.error("No INPUT files given")
        ctx.exit(2)

    if not name:
        name = inputs[0].name.split(".")[0]
        if not name:
            name = "input"

        logger.debug(f"No sample name found; using '{name}'")

    infile = tmp / f"{name}.fq.gz"
    if len(inputs) == 1 and inputs[0].is_file():
        if not inputs[0].suffix == ".gz":
            with open(inputs[0], "rb") as f_in:
                with gzip.open(infile, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
        else:
            shutil.copy2(inputs[0], infile)

    else:
        logger.info("Searching for fastq files...")
        fq_files = []
        for p in inputs:
            if p.is_file():
                fq_files.append(p)
            else:
                fq_files.extend(find_fastq_files(p, recursive))
        logger.info(f"Found {len(fq_files)} fastq files. Joining them...")
        concatenate_fastqs(fq_files, infile)

    # todo: download the reference genomes from https://www.nature.com/articles/srep45258#Sec25 and https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0214088
    # todo: remove contamination
    # todo: run mykrobe
    # todo: align decontam to h37rv
    # todo: bcftools mpileup/call
    # todo: filter vcf
    mykrobe = ExternalTool(
        tool="mykrobe",
        input=f"-i {infile}",
        output=f"-o {tmp}/{name}.mykrobe.json",
        params=f"predict --ont --force -A --debug --format json --min_proportion_expected_depth 0.20 --sample {name} "
               f"--species tb -t {threads} -m {mem_mb}MB --tmp {tmp} --skeleton_dir {tmp}"
    )

    subsampled_reads = f"{tmp}/{name}.subsampled.fastq.gz"
    rasusa = ExternalTool(
        tool="rasusa",
        input=f"-i {infile}",
        output=f"-o {subsampled_reads}",
        params=f"-c 150 -g 4411532 -s 88"
    )

    sam_file = f"{tmp}/{name}.subsampled.sam"
    minimap = ExternalTool(
        tool="minimap2",
        input=f"{H37RV_genome} {subsampled_reads}",
        output=f"-o {sam_file}",
        params=f"-a -L --sam-hit-only --secondary=no -x map-ont -t {threads}"
    )

    sorted_sam_file = f"{tmp}/{name}.subsampled.sorted.sam"
    samtools_sort = ExternalTool(
        tool="samtools",
        input=sam_file,
        output=f"-o {sorted_sam_file}",
        params=f"sort -@ {threads}"
    )

    pileup_file = f"{tmp}/{name}.subsampled.pileup.bcf"
    bcftools_mpileup = ExternalTool(
        tool="bcftools",
        input=sorted_sam_file,
        output=f"-o {pileup_file}",
        params=f"mpileup -x -O b -I -Q 13 -a 'INFO/SCR,FORMAT/SP,INFO/ADR,INFO/ADF' -h100 -M10000 -f {H37RV_genome} "
               f"--threads {threads}"
    )

    snps_file = f"{tmp}/{name}.subsampled.snps.bcf"
    bcftools_call = ExternalTool(
        tool="bcftools",
        input=pileup_file,
        output=f"-o {snps_file}",
        params=f"call --ploidy 1 -O b -V indels -m --threads {threads}"
    )

    tools_to_run = [mykrobe, rasusa, minimap, samtools_sort, bcftools_mpileup, bcftools_call]
    for tool in tools_to_run:
        try:
            tool.run()
        except subprocess.CalledProcessError as error:
            logger.error(f"Error calling {tool.command} (return code {error.returncode})")
            logger.error(f"Please check stdout log file: {tool.out_log}")
            logger.error(f"Please check stderr log file: {tool.err_log}")
            logger.error(f"Temporary files are preserved for debugging")
            logger.error("Exiting...")
            ctx.exit(1)

    if cleanup:
        logger.info("Cleaning up temporary files...")
        shutil.rmtree(tmp)

    logger.success("Done")


if __name__ == "__main__":
    main()
