import sys
from pathlib import Path

import click
from loguru import logger

from tbpore import __version__
from tbpore.cli import Mutex

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
@click.pass_context
def main(ctx: click.Context, verbose: bool, quiet: bool, outdir: Path):
    """Mycobacterium tuberculosis genomic analysis from Nanopore sequencing data"""
    log_lvl = "INFO"
    if verbose:
        log_lvl = "DEBUG"
    elif quiet:
        log_lvl = "ERROR"
    logger.remove()
    logger.add(sys.stderr, level=log_lvl, format=log_fmt)
    logger.info(f"Welcome to TBpore version {__version__}")


if __name__ == "__main__":
    main()
