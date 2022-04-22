# from https://raw.githubusercontent.com/mbhall88/head_to_head_pipeline/master/analysis/baseline_variants/scripts/consensus.py
# with a small modification to read gzipped references

import gzip
import logging
from collections import defaultdict
from enum import Enum
from pathlib import Path
from typing import Dict, List, NamedTuple, Optional, Set, TextIO, Union

import click
from cyvcf2 import VCF, Variant

N = "N"
PathLike = Union[str, Path]
DEFAULT_HET_CHOICES = ["none", "ref", "alt"]


class UnknownDefaultHet(Exception):
    pass


class Bed:
    """Positions in this class are 1-based - in contrast to BED positions being 0-based.
    This is to align it with VCF positions which are 1-based."""

    def __init__(self, file: Optional[PathLike] = None, zero_based: bool = True):
        self.file = Path(file) if file is not None else file
        self.zero_based = zero_based
        self.positions: Dict[str, Set[int]] = defaultdict(set)
        if self.file is not None:
            with self.file.open() as instream:
                for line in map(str.rstrip, instream):
                    # start is 0-based inclusive; end is 0-based non-inclusive
                    chrom, start, end = line.split()[:3]
                    start = int(start)
                    end = int(end)
                    self.positions[chrom].update(range(start, end))

    def __contains__(self, variant: Variant) -> bool:
        chrom = variant.CHROM
        pos = variant.POS if not self.zero_based else variant.POS - 1
        return chrom in self.positions and pos in self.positions[chrom]


class Genotype(NamedTuple):
    allele1: int
    allele2: int

    def is_null(self) -> bool:
        """Is the genotype null. i.e. ./."""
        return self.allele1 == -1 and self.allele2 == -1

    def is_hom(self) -> bool:
        """Is the genotype homozygous"""
        if self.is_null():
            return False
        if self.allele1 == -1 or self.allele2 == -1:
            return True
        return self.allele1 == self.allele2

    def is_het(self) -> bool:
        """Is the genotype heterozyhous"""
        return not self.is_null() and not self.is_hom()

    def is_hom_ref(self) -> bool:
        """Is genotype homozygous reference?"""
        return self.is_hom() and (self.allele1 == 0 or self.allele2 == 0)

    def is_hom_alt(self) -> bool:
        """Is genotype homozygous alternate?"""
        return self.is_hom() and (self.allele1 > 0 or self.allele2 > 0)

    def alt_index(self) -> Optional[int]:
        """If the genotype is homozygous alternate, returns the 0-based index of the
        alt allele in the alternate allele array.
        """
        if not self.is_hom_alt():
            return None
        return max(self.allele1, self.allele2) - 1

    @staticmethod
    def from_arr(arr: List[int]) -> "Genotype":
        alleles = [a for a in arr if type(a) is int]
        if len(alleles) < 2:
            alleles.append(-1)
        return Genotype(*alleles)


class Classification(Enum):
    Ref = "REF"
    Alt = "ALT"
    Null = "NULL"
    Het = "HET"

    def __str__(self):
        return self.value

    @staticmethod
    def from_variant(variant: Variant) -> "Classification":
        gt = Genotype.from_arr(variant.genotypes[0])
        if gt.is_null():
            classification = Classification.Null
        elif gt.is_het():
            classification = Classification.Het
        elif gt.is_hom_ref():
            classification = Classification.Ref
        elif gt.is_hom_alt():
            classification = Classification.Alt
        else:
            raise NotImplementedError(
                f"Can't determine genotype for variant \n{str(variant)}"
            )
        return classification


class Classifier:
    def __init__(
        self,
        mask: Optional[Bed] = None,
        ignore_filter: bool = False,
        ignore_mask: bool = False,
        ignore_null: bool = False,
        het_default: str = "none",
    ):
        self.ignore_filter = ignore_filter
        self.ignore_mask = ignore_mask
        self.ignore_null = ignore_null
        if het_default not in DEFAULT_HET_CHOICES:
            raise UnknownDefaultHet(
                f"{het_default} is not in the known het default choices: "
                f"{DEFAULT_HET_CHOICES}"
            )
        self.het_default = het_default

        if mask is None:
            mask = Bed()
        self.mask = mask

    def classify(self, variant: Variant) -> str:
        length = len(variant.REF)
        if self.ignore_mask and variant in self.mask:
            return N * length

        if self.ignore_filter and variant.FILTER is not None:
            return N * length

        classification = Classification.from_variant(variant)
        if classification is Classification.Null:
            return N * length if self.ignore_null else variant.REF
        if classification is Classification.Het:
            if self.het_default == "none":
                return N * length
            else:
                return_ref = self.het_default == "ref"
                if return_ref:
                    return variant.REF
                ref_in_called_alleles = 0 in variant.genotypes[0]
                if not ref_in_called_alleles:
                    logging.warning(
                        "Got a het with two different ALT calls. You have specified "
                        "ALT as the default base for hets. Using ALT with highest "
                        "index..."
                    )
                return variant.ALT[max(variant.genotypes[0]) - 1]

        if classification is Classification.Ref:
            return variant.REF
        if classification is Classification.Alt:
            gt = Genotype.from_arr(variant.genotypes[0])
            idx = gt.alt_index()
            assert idx is not None
            return variant.ALT[idx]
        raise NotImplementedError(f"Unable to classify variant: {str(variant)}")


def load_reference(path: str) -> Dict[str, List[str]]:
    index = dict()
    with gzip.open(path, "rt") as instream:
        seqid = None
        for line in map(str.rstrip, instream):
            if line.startswith(">"):
                seqid = line.split()[0][1:]
                if seqid in index:
                    raise ReferenceError(
                        f"Found duplicate sequence identifier in the reference fasta: "
                        f"{seqid}"
                    )
                index[seqid] = []
                continue
            if seqid is None:
                raise ReferenceError(
                    "Did not find header before sequence in reference fasta"
                )
            index[seqid].extend(list(line))

    return index


@click.command()
@click.help_option("--help", "-h")
@click.option(
    "-i",
    "--vcf",
    help="VCF file generate the consensus for.",
    type=click.Path(exists=True, dir_okay=False, allow_dash=True),
    default="-",
    show_default=True,
)
@click.option(
    "-f",
    "--ref",
    help="Reference fasta file.",
    type=click.Path(exists=True, dir_okay=False, allow_dash=False),
    required=True,
)
@click.option(
    "-m",
    "--mask",
    "bedfile",
    help="BED file containing positions to mask.",
    type=click.Path(exists=True, dir_okay=False),
)
@click.option(
    "-o",
    "--output",
    help="File path of consensus fasta to output.",
    type=click.File(mode="w"),
    default="-",
    show_default=True,
)
@click.option(
    "-H",
    "--het-default",
    help=(
        "Which heterozygous base should be used? If 'none', an 'N' will be used in the "
        "consensus sequence"
    ),
    type=click.Choice(DEFAULT_HET_CHOICES, case_sensitive=False),
    default=DEFAULT_HET_CHOICES[0],
    show_default=True,
)
@click.option(
    "-I",
    "--ignore",
    help=(
        "Which types of variants should be ignored (i.e. replaced with 'N')? If not "
        "ignoring filters/mask and the position is called ALT, then ALT will be used. "
        "For null/missing the REF base will be used. This option can be specified "
        "multiple times. For example, to ignore null and mask positions use: -I null "
        "-I mask"
    ),
    type=click.Choice(
        ("none", "all", "missing", "null", "filter", "mask"), case_sensitive=False
    ),
    default=["all"],
    show_default=True,
    multiple=True,
)
@click.option(
    "-s",
    "--sample-id",
    help=(
        "Sample identifier to use for the output fasta file. If not specified, the "
        "sample column name in the VCF will be used. If more than one contig is "
        "present, it will be combined with the sample ID (i.e. CONTIG|SAMPLE_ID)"
    ),
)
@click.option("-v", "--verbose", help="Turns on debug-level logging.", is_flag=True)
def main(
    vcf: str,
    ref: str,
    output: TextIO,
    verbose: bool,
    bedfile: str,
    het_default: str,
    ignore: tuple,
    sample_id: str,
):
    """This script will produce a consensus fasta file from a VCF file.
    The VCF must contain ONLY SNVs; at this stage it cannot handle any other variant
    type. It will also only use the first sample in the VCF currently.
    """
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s]: %(message)s", level=log_level
    )

    if "none" in ignore and len(ignore) > 1:
        raise click.BadOptionUsage(
            "ignore",
            f"'none' cannot be used in conjunction with other --ignore options and you "
            f"passed {ignore}",
        )

    if "mask" in ignore and bedfile is None:
        logging.warning(
            "You have asked to ignore masked positions, but no mask was provided"
        )

    if bedfile:
        logging.info("Loading mask...")
        mask = Bed(bedfile)
        logging.info(f"Loaded mask for {len(mask.positions)} chromosome(s).")
    else:
        logging.info("No mask given...")
        mask = Bed()

    logging.info("Loading reference fasta...")
    consensus = load_reference(ref)
    expected_consensus_seq_lens = {k: len(s) for k, s in consensus.items()}
    logging.info(f"Loaded {len(consensus)} sequence(s) from reference fasta")

    vcf_reader = VCF(vcf)

    if not sample_id:
        sample_id = vcf_reader.samples[0]

    multiple_contigs = len(vcf_reader.seqnames) > 1
    if multiple_contigs:
        logging.info(
            f"Sequence identifiers will be CONTIG|{sample_id} where CONTIG is the name "
            f"of the relevant contig"
        )
    else:
        logging.info(f"Sequence identifier will be {sample_id}")

    ignore_filter = any(s in ignore for s in ("filter", "all"))
    ignore_null = any(s in ignore for s in ("null", "all"))
    ignore_missing = any(s in ignore for s in ("missing", "all"))
    ignore_mask = any(s in ignore for s in ("mask", "all"))
    classifier = Classifier(
        mask=mask,
        ignore_filter=ignore_filter,
        ignore_mask=ignore_mask,
        ignore_null=ignore_null,
        het_default=het_default,
    )
    warned_chroms = set()
    positions_in_vcf = defaultdict(set)
    logging.info("Generating consensus from variants...")
    for variant in vcf_reader:
        if variant.CHROM not in consensus and variant.CHROM not in warned_chroms:
            logging.warning(
                f"{variant.CHROM} is not in the reference fasta. Skipping..."
            )
            warned_chroms.add(variant.CHROM)
            continue

        bases = classifier.classify(variant)
        start_on_consensus = variant.POS - 1
        end_on_consensus = start_on_consensus + len(bases)
        consensus[variant.CHROM][start_on_consensus:end_on_consensus] = list(bases)
        for i in range(start_on_consensus, end_on_consensus):
            positions_in_vcf[variant.CHROM].add(i)
    logging.info("Consensus from variants complete")

    if ignore_missing:
        logging.info("Updating consensus with ignored missing positions...")
        for chrom in consensus:
            all_positions = set(range(len(consensus[chrom])))
            missing_positions = all_positions - positions_in_vcf.get(chrom, set())
            for pos in missing_positions:
                consensus[chrom][pos] = N
        logging.info("Finished updating missing positions")

    logging.info("Writing output consensus fasta...")
    for chrom in consensus:
        if not multiple_contigs:
            header = f">{sample_id}"
        else:
            header = f">{chrom}|{sample_id}"
        print(header, file=output)
        consensus_seq = "".join(consensus[chrom])
        expected_length = expected_consensus_seq_lens[chrom]
        assert (
            len(consensus_seq) == expected_length
        ), f"{chrom}'s consensus sequence ({len(consensus_seq)}) is not the expected length {expected_length}"
        print(consensus_seq, file=output)

    logging.info("Done")


if __name__ == "__main__":
    main()
