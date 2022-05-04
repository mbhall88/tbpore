# From https://github.com/mbhall88/head_to_head_pipeline/blob/master/analysis/baseline_variants/scripts/apply_filters.py

import logging
from collections import Counter
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import List, NamedTuple, Optional, Tuple

import click
import numpy as np
from cyvcf2 import VCF, Variant, Writer

HIST_BINS = 40


def validate_fraction(ctx, param, value):
    f = float(value)
    if f > 1.0 or f < 0.0:
        raise click.BadParameter("fractions must be between 0.0 and 1.0")
    else:
        return f


class Tags(Enum):
    Depth = "DP"
    LowDepth = "ld"
    LowFed = "lfed"
    HighDepth = "hd"
    LowQual = "lq"
    StrandBias = "sb"
    StrandDepth = "DP4"
    Pass = "PASS"
    BaseQualBias = "BQB"
    LowBaseQualBias = "lbqb"
    MapQualBias = "MQB"
    LowMapQualBias = "lmqb"
    ReadPosBias = "RPB"
    ReadPosBiasZ = "RPBZ"
    LowReadPosBias = "lrpb"
    LowReadPosBiasZ = "lrpbz"
    HighReadPosBiasZ = "hrpbz"
    SegregationBias = "SGB"
    HighSegBias = "hsgb"
    SoftClipBiasZ = "SCBZ"
    HighSoftClipBiasZ = "hscbz"
    VariantDistanceBias = "VDB"
    LowVarDistBias = "lvdb"
    LowSupport = "frs"
    LowMapQual = "lmq"
    MapQual = "MQ"

    def __str__(self) -> str:
        return str(self.value)


@dataclass
class FilterStatus:
    low_depth: bool = False
    low_fed: bool = False
    high_depth: bool = False
    low_qual: bool = False
    strand_bias: bool = False
    low_bqb: bool = False
    low_mqb: bool = False
    low_rpb: bool = False
    low_rpbz: bool = False
    high_rpbz: bool = False
    high_scbz: bool = False
    low_mapq: bool = False
    high_sgb: bool = False
    low_vdb: bool = False
    low_support: bool = False
    delim: str = ";"

    def __str__(self) -> str:
        status = []
        if self.low_depth:
            status.append(str(Tags.LowDepth))
        if self.low_fed:
            status.append(str(Tags.LowFed))
        if self.high_depth:
            status.append(str(Tags.HighDepth))
        if self.low_qual:
            status.append(str(Tags.LowQual))
        if self.strand_bias:
            status.append(str(Tags.StrandBias))
        if self.low_bqb:
            status.append(str(Tags.LowBaseQualBias))
        if self.low_mqb:
            status.append(str(Tags.LowMapQualBias))
        if self.low_mapq:
            status.append(str(Tags.LowMapQual))
        if self.low_rpb:
            status.append(str(Tags.LowReadPosBias))
        if self.low_rpbz:
            status.append(str(Tags.LowReadPosBiasZ))
        if self.high_rpbz:
            status.append(str(Tags.HighReadPosBiasZ))
        if self.high_scbz:
            status.append(str(Tags.HighSoftClipBiasZ))
        if self.high_sgb:
            status.append(str(Tags.HighSegBias))
        if self.low_vdb:
            status.append(str(Tags.LowVarDistBias))
        if self.low_support:
            status.append(str(Tags.LowSupport))

        return self.delim.join(status) if status else str(Tags.Pass)


@dataclass
class StrandDepths:
    ref_forward: int = 0
    ref_reverse: int = 0
    alt_forward: int = 0
    alt_reverse: int = 0

    @staticmethod
    def _ratio(depths: Tuple[int, int]) -> float:
        try:
            return min(depths) / sum(depths)
        except ZeroDivisionError:
            return 1.0

    @property
    def ref_depths(self) -> Tuple[int, int]:
        return self.ref_forward, self.ref_reverse

    @property
    def alt_depths(self) -> Tuple[int, int]:
        return self.alt_forward, self.alt_reverse

    @property
    def ref_ratio(self) -> float:
        return self._ratio(self.ref_depths)

    @property
    def alt_ratio(self) -> float:
        return self._ratio(self.alt_depths)

    def to_tuple(self) -> Tuple[int, int, int, int]:
        return self.ref_forward, self.ref_reverse, self.alt_forward, self.alt_reverse


class DepthTagError(Exception):
    pass


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

    def allele_index(self) -> Optional[int]:
        """The index of the called allele"""
        if self.is_hom_ref() or self.is_null():
            return 0
        elif self.is_hom_alt():
            return self.alt_index() + 1
        else:
            raise NotImplementedError(f"Het Genotype is unexpected: {self}")

    @staticmethod
    def from_arr(arr: List[int]) -> "Genotype":
        alleles = [a for a in arr if type(a) is int]
        if len(alleles) < 2:
            alleles.append(-1)
        return Genotype(*alleles)


class Filter:
    def __init__(
        self,
        expected_depth: int = 0,
        min_depth: int = 0,
        min_fed: float = 0,
        max_depth: int = 0,
        min_strand_bias: int = 0,
        min_qual: float = 0,
        min_bqb: float = 0,
        min_mqb: float = 0,
        min_rpb: float = 0,
        max_sgb: float = 0,
        min_vdb: float = 0,
        min_frs: float = 0,
        min_mq: int = 0,
        min_rpbz: Optional[float] = None,
        max_rpbz: Optional[float] = None,
        max_scbz: Optional[float] = None,
    ):
        self.expected_depth = expected_depth
        self.min_depth = min_depth
        self.min_fed = min_fed
        self.max_depth = max_depth
        self.min_strand_bias = min_strand_bias / 100
        self.min_qual = min_qual
        self.min_bqb = min_bqb
        self.min_mqb = min_mqb
        self.min_mq = min_mq
        self.min_rpb = min_rpb
        self.max_sgb = max_sgb
        self.min_vdb = min_vdb
        self.min_frs = min_frs
        if min_rpbz is None:
            min_rpbz = -float("inf")
        self.min_rpbz = min_rpbz
        if max_rpbz is None:
            max_rpbz = float("inf")
        self.max_rpbz = max_rpbz
        if max_scbz is None:
            max_scbz or float("inf")
        self.max_scbz = max_scbz

        if self.min_depth and self.max_depth and self.min_depth > self.max_depth:
            raise ValueError(
                f"Minimum depth is more than maximum depth: "
                f"{self.min_depth:.1f} > {self.max_depth:.1f}"
            )

    def _is_low_depth(self, variant: Variant) -> bool:
        variant_depth = get_depth(variant)
        return variant_depth < self.min_depth if self.min_depth else False

    def _is_low_fed(self, variant: Variant) -> bool:
        hq_depth = get_hq_depth(variant)
        fed = hq_depth / self.expected_depth
        return fed < self.min_fed if self.min_fed > 0 else False

    def _is_high_depth(self, variant: Variant) -> bool:
        variant_depth = get_depth(variant)
        return variant_depth > self.max_depth if self.max_depth else False

    def _is_low_qual(self, variant: Variant) -> bool:
        return (variant.QUAL or 0) < self.min_qual

    def _is_low_mapq(self, variant: Variant) -> bool:
        variant_mapq = get_mapq(variant)
        return variant_mapq < self.min_mq if self.min_mq else False

    def _is_low_bqb(self, variant: Variant) -> bool:
        bqb = variant.INFO.get(str(Tags.BaseQualBias))
        return bqb is not None and bqb < self.min_bqb

    def _is_low_mqb(self, variant: Variant) -> bool:
        mqb = variant.INFO.get(str(Tags.MapQualBias))
        return mqb is not None and mqb < self.min_mqb

    def _is_low_rpb(self, variant: Variant) -> bool:
        rpb = variant.INFO.get(str(Tags.ReadPosBias))
        return rpb is not None and rpb < self.min_rpb

    def _is_high_sgb(self, variant: Variant) -> bool:
        sgb = variant.INFO.get(str(Tags.SegregationBias))
        return sgb is not None and sgb > self.max_sgb

    def _is_low_vdb(self, variant: Variant) -> bool:
        vdb = variant.INFO.get(str(Tags.VariantDistanceBias))
        return vdb is not None and vdb < self.min_vdb

    def _is_low_support(self, variant: Variant) -> bool:
        frs = fraction_read_support(variant)
        return frs < self.min_frs if self.min_frs else False

    def _is_low_rpbz(self, variant: Variant) -> bool:
        rpbz = variant.INFO.get(str(Tags.ReadPosBiasZ))
        return rpbz is not None and rpbz < self.min_rpbz

    def _is_high_rpbz(self, variant: Variant) -> bool:
        rpbz = variant.INFO.get(str(Tags.ReadPosBiasZ))
        return rpbz is not None and rpbz > self.max_rpbz

    def _is_high_scbz(self, variant: Variant) -> bool:
        scbz = variant.INFO.get(str(Tags.SoftClipBiasZ))
        return scbz is not None and scbz > self.max_scbz

    def filter_status(self, variant: Variant) -> str:
        status = FilterStatus()
        if self.min_depth or self.max_depth:
            status.low_depth = self._is_low_depth(variant)
            status.high_depth = self._is_high_depth(variant)

        if self.min_fed > 0:
            status.low_fed = self._is_low_fed(variant)

        if self.min_mq:
            status.low_mapq = self._is_low_mapq(variant)

        if self.min_qual:
            status.low_qual = self._is_low_qual(variant)

        if self.min_strand_bias:
            strand_depths = get_strand_depths(variant)
            assert strand_depths is not None, (
                f"Strand bias filter should be turned off if no {Tags.StrandDepth} "
                f"tag is present."
            )

            gt = Genotype.from_arr(variant.genotypes[0])
            if gt.is_hom_alt():
                ratio = strand_depths.alt_ratio
            elif gt.is_hom_ref():
                ratio = strand_depths.ref_ratio
            elif gt.is_het():
                ratio = min(strand_depths.ref_ratio, strand_depths.alt_ratio)
            elif gt.is_null():
                ratio = float("inf")
            else:
                raise NotImplementedError(
                    f"Don't know how to interpret genotype {gt} for variant at "
                    f"POS {variant.POS}"
                )
            status.strand_bias = ratio < self.min_strand_bias

        if self.min_bqb:
            status.low_bqb = self._is_low_bqb(variant)

        if self.min_mqb:
            status.low_mqb = self._is_low_mqb(variant)

        if self.min_rpb:
            status.low_rpb = self._is_low_rpb(variant)

        if self.min_rpbz:
            status.low_rpbz = self._is_low_rpbz(variant)

        if self.max_rpbz:
            status.high_rpbz = self._is_high_rpbz(variant)

        if self.max_scbz:
            status.high_scbz = self._is_high_scbz(variant)

        if self.max_sgb != 0:
            status.high_sgb = self._is_high_sgb(variant)

        if self.min_vdb:
            status.low_vdb = self._is_low_vdb(variant)

        if self.min_frs:
            status.low_support = self._is_low_support(variant)

        return str(status)

    def add_filters_to_header(self, vcf: VCF):
        if self.min_depth > 0:
            header = {
                "ID": str(Tags.LowDepth),
                "Description": (
                    f"Depth ({Tags.Depth}) less than {self.min_depth} - i.e., {Tags.Depth}<{self.min_depth:.1f}"
                ),
            }
            vcf.add_filter_to_header(header)
            logging.debug(f"Header for min. depth: {header}")

        if self.min_fed > 0:
            header = {
                "ID": str(Tags.LowFed),
                "Description": (
                    f"High-quality depth of the called allele as a fraction of the expected (median; {self.expected_depth}) is "
                    f"less than {self.min_fed}"
                ),
            }
            vcf.add_filter_to_header(header)
            logging.debug(f"Header for min. FED: {header}")

        if self.min_mq > 0:
            header = {
                "ID": str(Tags.LowMapQual),
                "Description": (
                    f"Mapping quality ({Tags.MapQual.value}) less than {self.min_mq} - i.e., {Tags.MapQual.value}<{self.min_mq:.0f}"
                ),
            }
            vcf.add_filter_to_header(header)
            logging.debug(f"Header for min. depth: {header}")

        if self.max_depth > 0:
            header = {
                "ID": str(Tags.HighDepth),
                "Description": (
                    f"Depth ({Tags.Depth}) more than {self.max_depth} - i.e., {Tags.Depth}>{self.max_depth:.1f}"
                ),
            }
            vcf.add_filter_to_header(header)
            logging.debug(f"Header for max. depth: {header}")

        if self.min_qual > 0:
            header = {
                "ID": str(Tags.LowQual),
                "Description": f"QUAL less than {self.min_qual}",
            }
            vcf.add_filter_to_header(header)
            logging.debug(f"Header for min. QUAL: {header}")

        if self.min_strand_bias > 0:
            header = {
                "ID": str(Tags.StrandBias),
                "Description": (
                    f"A strand on the called allele has less than  "
                    f"{self.min_strand_bias:.2%} of the high-quality depth for that "
                    f"allele. This is judged on the {Tags.StrandDepth} tag."
                ),
            }
            vcf.add_filter_to_header(header)
            logging.debug(f"Header for strand bias: {header}")

        if self.min_frs > 0:
            header = {
                "ID": str(Tags.LowSupport),
                "Description": f"Fraction of read support on called allele is less than {self.min_frs}",
            }
            vcf.add_filter_to_header(header)
            logging.debug(f"Header for min. FRS: {header}")

        if self.min_bqb > 0:
            header = {
                "ID": str(Tags.LowBaseQualBias),
                "Description": (
                    f"Base Quality Bias ({Tags.BaseQualBias}) is less than "
                    f"{self.min_bqb}."
                ),
            }
            vcf.add_filter_to_header(header)
            logging.debug(f"Header for min. base quality bias: {header}")

        if self.min_mqb > 0:
            header = {
                "ID": str(Tags.LowMapQualBias),
                "Description": (
                    f"Mapping Quality Bias ({Tags.MapQualBias}) is less than "
                    f"{self.min_mqb}."
                ),
            }
            vcf.add_filter_to_header(header)
            logging.debug(f"Header for min. mapping quality bias: {header}")

        if self.min_rpb > 0:
            header = {
                "ID": str(Tags.LowReadPosBias),
                "Description": (
                    f"Read Position Bias ({Tags.ReadPosBias}) is less than "
                    f"{self.min_rpb}."
                ),
            }
            vcf.add_filter_to_header(header)
            logging.debug(f"Header for min. read position bias: {header}")

        if self.min_rpbz is not None:
            header = {
                "ID": str(Tags.LowReadPosBiasZ),
                "Description": (
                    f"Read Position Bias z-test score ({Tags.ReadPosBiasZ}) is less than "
                    f"{self.min_rpbz}."
                ),
            }
            vcf.add_filter_to_header(header)
            logging.debug(f"Header for min. read position bias z-test: {header}")

        if self.max_rpbz is not None:
            header = {
                "ID": str(Tags.HighReadPosBiasZ),
                "Description": (
                    f"Read Position Bias z-test score ({Tags.ReadPosBiasZ}) is more than "
                    f"{self.max_rpbz}."
                ),
            }
            vcf.add_filter_to_header(header)
            logging.debug(f"Header for max. read position bias z-test: {header}")

        if self.max_scbz is not None:
            header = {
                "ID": str(Tags.HighSoftClipBiasZ),
                "Description": (
                    f"Soft-Clip Length Bias z-test score ({Tags.SoftClipBiasZ}) is more than "
                    f"{self.max_scbz}."
                ),
            }
            vcf.add_filter_to_header(header)
            logging.debug(f"Header for max. soft-clip length bias z-test: {header}")

        if self.max_sgb != 0:
            header = {
                "ID": str(Tags.HighSegBias),
                "Description": (
                    f"Segregation-based metric ({Tags.SegregationBias}) is greater "
                    f"than {self.max_sgb}."
                ),
            }
            vcf.add_filter_to_header(header)
            logging.debug(f"Header for max. segregation bias: {header}")

        if self.min_vdb > 0:
            header = {
                "ID": str(Tags.LowVarDistBias),
                "Description": (
                    f"Variant distance bias ({Tags.VariantDistanceBias}) is less "
                    f"than {self.min_vdb}."
                ),
            }
            vcf.add_filter_to_header(header)
            logging.debug(f"Header for min. variant distance bias: {header}")


def fraction_read_support(variant: Variant, sample_idx: int = 0) -> float:
    strand_depths = get_strand_depths(variant)
    if strand_depths is None:
        return 1.0
    total_covg = sum(strand_depths.to_tuple())
    called_idx = Genotype.from_arr(variant.genotypes[0]).allele_index()
    if called_idx == 0:  # ref
        called_covg = sum(strand_depths.ref_depths)
    else:
        called_covg = sum(strand_depths.alt_depths)
    try:
        return called_covg / total_covg
    except ZeroDivisionError:
        return 1.0


def get_depth(variant: Variant, default: int = 0) -> int:
    return variant.INFO.get(str(Tags.Depth), default)


def get_hq_depth(variant: Variant, default: int = 0) -> int:
    strand_depths = get_strand_depths(variant)
    if strand_depths is None:
        return default
    called_idx = Genotype.from_arr(variant.genotypes[0]).allele_index()
    if called_idx == 0:  # ref
        called_covg = sum(strand_depths.ref_depths)
    else:
        called_covg = sum(strand_depths.alt_depths)
    return called_covg


def get_mapq(variant: Variant, default: int = 0) -> int:
    mq = variant.INFO.get(str(Tags.MapQual), default)
    if mq is None:  # even though we gave a default it can still return none
        mq = default
    return mq


def get_strand_depths(
    variant: Variant, default: Optional[StrandDepths] = None
) -> Optional[StrandDepths]:
    strand_depths = variant.INFO.get(str(Tags.StrandDepth), None)
    return StrandDepths(*strand_depths) if strand_depths is not None else default


@click.command()
@click.help_option("--help", "-h")
@click.option(
    "-i",
    "--in-vcf",
    help="VCF file to apply filters to.",
    type=click.Path(exists=True, dir_okay=False, allow_dash=True),
    default="-",
    show_default=True,
)
@click.option(
    "-o",
    "--out-vcf",
    help="New VCF with filters.",
    type=click.Path(exists=False, dir_okay=False, allow_dash=True),
    default="-",
    show_default=True,
)
@click.option(
    "-d",
    "--min-depth",
    help=(
        "Minimum read depth. This filter "
        f"has ID: {Tags.LowDepth.value}. Set to 0 to disable"
    ),
    default=0,
    show_default=True,
)
@click.option(
    "-x",
    "--min-fed",
    help=f"Minimum fraction of expected depth. Where expected depth is the median. "
    f"This filter has ID: {Tags.LowFed.value}. Set to 0 to disable",
    default=0.0,
    show_default=True,
    type=float,
)
@click.option(
    "-D",
    "--max-depth",
    help=(
        "Maximum depth. This filter "
        f"has ID: {Tags.HighDepth.value}. Set to 0 to disable"
    ),
    default=0,
    show_default=True,
)
@click.option(
    "-s",
    "--min-strand-bias",
    help=(
        "Filter a variant if either strand has less than INT% of the (high-quality) "
        f"depth on the called allele ({Tags.StrandDepth.value}). This filter has ID: "
        f"{Tags.StrandBias}. Set to 0 to disable"
    ),
    type=click.IntRange(0, 50),
    metavar="INT",
    default=25,
    show_default=True,
)
@click.option(
    "-q",
    "--min-qual",
    help=(
        f"Filter a variant if QUAL is less than INT. This filter has ID: "
        f"{Tags.LowQual.value}. Set to 0 to disable"
    ),
    default=0.0,
    show_default=True,
)
@click.option(
    "-K",
    "--min-frs",
    help=(
        f"Minimum fraction of reads supporting the called allele. "
        f"This filter has ID: {Tags.LowSupport.value}. Set to 0 to disable"
    ),
    default=0.0,
    type=float,
    show_default=True,
    callback=validate_fraction,
)
@click.option(
    "-M",
    "--min-mq",
    help=(
        f"Minimum mapping quality ({Tags.MapQual.value}) score. "
        f"This filter has ID: {Tags.LowMapQual.value}. Set to 0 to disable"
    ),
    default=0,
    type=int,
    show_default=True,
)
@click.option(
    "-b",
    "--min-bqb",
    help=(
        f"Filter a variant if base quality bias ({Tags.BaseQualBias.value}) is less "
        f"than FLOAT. This filter has ID: {Tags.LowBaseQualBias.value}. Set to 0 to "
        f"disable"
    ),
    default=0.0,
    show_default=True,
)
@click.option(
    "-m",
    "--min-mqb",
    help=(
        f"Filter a variant if mapping quality bias ({Tags.MapQualBias.value}) is "
        f"less than FLOAT. This filter has ID: {Tags.LowMapQualBias.value}. Set to 0 to "
        f"disable"
    ),
    default=0.0,
    show_default=True,
)
@click.option(
    "-r",
    "--min-rpb",
    help=(
        f"Filter a variant if read position bias ({Tags.ReadPosBias.value}) is "
        f"less than FLOAT. This filter has ID: {Tags.LowReadPosBias.value}. Set to 0 to "
        f"disable"
    ),
    default=0.0,
    show_default=True,
)
@click.option(
    "-w",
    "--min-rpbz",
    help=(
        f"Filter a variant if read position bias z-test score ({Tags.ReadPosBiasZ.value}) is "
        f"less than FLOAT. This filter has ID: {Tags.LowReadPosBiasZ.value}."
    ),
    default=None,
    type=float,
)
@click.option(
    "-W",
    "--max-rpbz",
    help=(
        f"Filter a variant if read position bias z-test score ({Tags.ReadPosBiasZ.value}) is "
        f"more than FLOAT. This filter has ID: {Tags.HighReadPosBiasZ.value}."
    ),
    default=None,
    type=float,
)
@click.option(
    "-C",
    "--max-scbz",
    help=(
        f"Filter a variant if soft-clip length bias z-test score ({Tags.SoftClipBiasZ.value}) is "
        f"more than FLOAT. This filter has ID: {Tags.HighSoftClipBiasZ.value}."
    ),
    default=None,
    type=float,
)
@click.option(
    "-G",
    "--max-sgb",
    help=(
        f"Filter a variant if segregation bias ({Tags.SegregationBias.value}) is "
        f"greater than FLOAT. This filter has ID: {Tags.HighSegBias.value}. Set to 0 to "
        f"disable"
    ),
    default=0.0,
    show_default=True,
)
@click.option(
    "-V",
    "--min-vdb",
    help=(
        f"Filter a variant if variant distance bias ({Tags.VariantDistanceBias.value}) "
        f"is less than FLOAT. This filter has ID: {Tags.LowVarDistBias.value}. Set to 0 "
        f"to disable"
    ),
    default=0.0,
    show_default=True,
)
@click.option(
    "--overwrite/--no-overwrite",
    "-f/-F",
    default=True,
    show_default=True,
    help="Overwrite existing information in FILTER field.",
)
@click.option(
    "-p/-P",
    "--hist/--no-hist",
    default=True,
    show_default=True,
    help="Print histograms of depth and QUAL. Histograms will be printed to stderr.",
)
@click.option("-v", "--verbose", help="Turns on debug-level logging.", is_flag=True)
def main(
    in_vcf: str,
    out_vcf: str,
    overwrite: bool,
    verbose: bool,
    min_qual: float,
    min_depth: int,
    min_fed: float,
    max_depth: int,
    min_strand_bias: int,
    min_bqb: float,
    min_mqb: float,
    min_rpb: float,
    min_rpbz: Optional[float],
    max_rpbz: Optional[float],
    max_scbz: Optional[float],
    max_sgb: float,
    min_vdb: float,
    hist: bool,
    min_frs: float,
    min_mq: int,
):
    """Apply the following filters to a VCF:\n
    - Minimum proportion of the expected (median) depth\n
    - Maximum proportion of the expected (median) depth\n
    - Minimum QUAL threshold\n
    - Minimum Strand bias percentage
    """
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s]: %(message)s", level=log_level
    )

    vcf_reader = VCF(in_vcf)
    if not vcf_reader.contains(Tags.Depth.value):
        raise DepthTagError(f"Depth tag {Tags.Depth} not found in header")

    if (not vcf_reader.contains(str(Tags.StrandDepth))) and min_strand_bias:
        logging.warning(
            f"Strand depth tag {Tags.StrandDepth} not found in header. "
            f"Turning off strand bias filter..."
        )
        min_strand_bias = 0

    logging.info("Calculating expected (median) depth...")
    depths = []
    quals = []
    for v in vcf_reader:
        depths.append(get_depth(v))
        quals.append(v.QUAL or 0)

    median_depth = np.median(depths)
    logging.info(f"Expected depth: {median_depth}")

    if hist:
        import histoprint

        tick_format = "% .1f"
        logging.info("Depth histogram:")
        histoprint.print_hist(
            np.histogram(depths, bins=HIST_BINS),
            title="Depth histogram",
            summary=True,
            tick_format=tick_format,
            file=click.get_text_stream("stderr"),
        )

        logging.info("QUAL histogram")
        histoprint.print_hist(
            np.histogram(quals, bins=HIST_BINS),
            title="QUAL histogram",
            summary=True,
            tick_format=tick_format,
            file=click.get_text_stream("stderr"),
        )

    assessor = Filter(
        expected_depth=int(median_depth),
        min_qual=min_qual,
        min_depth=min_depth,
        min_fed=min_fed,
        max_depth=max_depth,
        min_strand_bias=min_strand_bias,
        min_bqb=min_bqb,
        min_mqb=min_mqb,
        min_rpb=min_rpb,
        max_sgb=max_sgb,
        min_vdb=min_vdb,
        min_frs=min_frs,
        min_mq=min_mq,
        min_rpbz=min_rpbz,
        max_rpbz=max_rpbz,
        max_scbz=max_scbz,
    )

    vcf_reader = VCF(in_vcf)
    assessor.add_filters_to_header(vcf_reader)

    if not Path(out_vcf).parent.exists():
        Path(out_vcf).parent.mkdir(exist_ok=True, parents=True)

    vcf_writer = Writer(out_vcf, tmpl=vcf_reader)

    stats = Counter()
    logging.info("Filtering variants...")
    for variant in vcf_reader:
        filter_status = assessor.filter_status(variant)

        if (
            (not overwrite)
            and variant.FILTER is not None
            and filter_status != str(Tags.Pass)
        ):
            current_filter = variant.FILTER.rstrip(";")
            variant.FILTER = f"{current_filter};{filter_status}"
        else:
            variant.FILTER = filter_status

        vcf_writer.write_record(variant)

        stats.update(filter_status.split(";"))

    vcf_reader.close()
    vcf_writer.close()

    logging.info("FILTER STATISTICS")
    logging.info("=================")
    for filt, count in stats.items():
        logging.info(f"Filter: {filt}\tCount: {count}")

    logging.info("Done!")


if __name__ == "__main__":
    main()
