import click

from tbpore import __version__


@click.command()
@click.help_option("--help", "-h")
@click.version_option(__version__, "--version", "-V")
@click.pass_context
def main(ctx: click.Context):
    """Mycobacterium tuberculosis genomic analysis from Nanopore sequencing data"""
    pass


if __name__ == "__main__":
    main()
