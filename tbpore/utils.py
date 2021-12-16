import fileinput
import glob
import gzip
import re
from pathlib import Path
from typing import IO, List, Sequence, Union

FASTQ_REGEX = re.compile(r"\.f(ast)?q(\.gz)?$")
PathLike = Union[str, Path]


def is_fastq(p: PathLike) -> bool:
    return Path(p).is_file() and bool(FASTQ_REGEX.search(str(p)))


def find_fastq_files(src: Path, recursive: bool = False) -> List[Path]:
    return [p for p in glob.iglob(f"{src}/**", recursive=recursive) if is_fastq(p)]


def which_open(fname: PathLike, mode: str) -> IO:
    if str(fname).endswith(".gz"):
        return gzip.open(fname, mode=mode)
    return open(fname, mode=mode)


def concatenate_fastqs(files: Sequence[PathLike], dest: PathLike):
    with gzip.open(dest, compresslevel=6, mode="wb") as fout, fileinput.input(
        files, openhook=which_open, mode="rb"
    ) as fin:
        for line in fin:
            if line.rstrip():  # skip empty lines
                fout.write(line)
