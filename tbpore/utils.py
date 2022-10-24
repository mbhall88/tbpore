import fileinput
import glob
import gzip
import re
from pathlib import Path
from typing import IO, Any, Dict, Set, Union

FASTQ_REGEX = re.compile(r"\.f(ast)?q(\.gz)?$")
PathLike = Union[str, Path]


def is_fastq(p: PathLike) -> bool:
    return Path(p).is_file() and bool(FASTQ_REGEX.search(str(p)))


def find_fastq_files(src: Path, recursive: bool = False) -> Set[str]:
    return set(p for p in glob.iglob(f"{src}/**", recursive=recursive) if is_fastq(p))


def which_open(fname: PathLike, mode: str) -> IO:
    if str(fname).endswith(".gz"):
        return gzip.open(fname, mode=mode)
    return open(fname, mode=mode)


def concatenate_fastqs(files: Set[PathLike], dest: PathLike):
    with gzip.open(dest, compresslevel=6, mode="wb") as fout, fileinput.input(
        files, openhook=which_open, mode="rb"
    ) as fin:
        for line in fin:
            if line.rstrip():  # skip empty lines
                fout.write(line)


def parse_verbose_filter_params(filters_dict: Dict[Any, Any]) -> str:
    opts = [
        ("d", "min_depth"),
        ("D", "max_depth"),
        ("q", "min_qual"),
        ("s", "min_strand_bias"),
        ("b", "min_bqb"),
        ("m", "min_mqb"),
        ("r", "min_rpb"),
        ("V", "min_vdb"),
        ("G", "max_sgb"),
        ("K", "min_frs"),
        ("w", "min_rpbz"),
        ("W", "max_rpbz"),
        ("C", "max_scbz"),
        ("M", "min_mq"),
        ("x", "min_fed"),
    ]
    flags = []
    for op, key in opts:
        if key in filters_dict:
            flags.append(f"-{op} {filters_dict[key]}")
    return " ".join(flags)


def fastq_prefix(path: Union[str, Path]) -> str:
    fname = Path(path).name
    return FASTQ_REGEX.sub("", fname, count=1)
