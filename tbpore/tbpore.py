import shutil
import sys
from pathlib import Path
from typing import Tuple

import click
from loguru import logger

from tbpore import TMP_NAME, __version__
from tbpore.cli import Mutex
from tbpore.utils import concatenate_fastqs, find_fastq_files

log_fmt = (
    "[<green>{time:YYYY-MM-DD HH:mm:ss}</green>] <level>{level: <8}</level> | "
    "<level>{message}</level>"
)


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

    outdir.mkdir(exist_ok=True)
    if tmp is None:
        tmp = outdir / TMP_NAME
    tmp.mkdir(exist_ok=True)

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

    if cleanup:
        logger.info("Cleaning up temporary files...")
        shutil.rmtree(tmp)

    logger.success("Done")


if __name__ == "__main__":
    main()
