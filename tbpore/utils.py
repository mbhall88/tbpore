import fileinput
import glob
import gzip
import hashlib
import re
import shutil
import urllib.request
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


def download_file(url: str, filename: Path):
    urllib.request.urlretrieve(url=url, filename=filename)


def validate_sha256(file_path: Path, expected_hash: str) -> bool:
    hash_obj = hashlib.sha256()

    # process the file in (byte) chunks
    chunk_size = 1_024 * 10_000
    with open(file_path, "rb") as f:
        while True:
            chunk = f.read(chunk_size)
            if not chunk:
                break
            hash_obj.update(chunk)

    # get the hexadecimal representation of the sha256 hash
    actual_hash = hash_obj.hexdigest()

    return actual_hash == expected_hash


def decompress_file(
    compressed_file: Path, decompressed_file: Path, remove_compressed: bool = False
):
    """Decompress a gzip-compressed file
    `remove_compressed` indicates whether to delete the compressed file after
    successful decompression
    """
    with gzip.open(compressed_file, "rb") as f_in:
        with open(decompressed_file, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    if remove_compressed:
        compressed_file.unlink(missing_ok=True)
