import gzip
import tempfile
from pathlib import Path
from unittest import mock

from tbpore.utils import (
    decompress_file,
    download_file,
    fastq_prefix,
    find_fastq_files,
    is_fastq,
    validate_sha256,
)


class TestIsFastq:
    @mock.patch("tbpore.utils.Path.is_file")
    def test_is_text_file(self, mock_is_file):
        mock_is_file.return_value = True
        p = "foo.txt"
        assert not is_fastq(p)

    @mock.patch("tbpore.utils.Path.is_file")
    def test_shortform(self, mock_is_file):
        mock_is_file.return_value = True
        p = "dir/a.fq"
        assert is_fastq(p)

    @mock.patch("tbpore.utils.Path.is_file")
    def test_longform(self, mock_is_file):
        mock_is_file.return_value = True
        p = "a.fastq"
        assert is_fastq(p)

    @mock.patch("tbpore.utils.Path.is_file")
    def test_short_compressed(self, mock_is_file):
        mock_is_file.return_value = True
        p = "tmp/foo.fq.gz"
        assert is_fastq(p)

    @mock.patch("tbpore.utils.Path.is_file")
    def test_long_compressed(self, mock_is_file):
        mock_is_file.return_value = True
        p = "bar.fastq.gz"
        assert is_fastq(p)

    @mock.patch("tbpore.utils.Path.is_file")
    def test_typo(self, mock_is_file):
        mock_is_file.return_value = True
        p = "bar.fatq.gz"
        assert not is_fastq(p)

    @mock.patch("tbpore.utils.Path.is_file")
    def test_different_compression_ext(self, mock_is_file):
        p = "bar.fastq.xz"
        mock_is_file.return_value = True
        assert not is_fastq(p)


@mock.patch("tbpore.utils.Path.is_file")
@mock.patch("tbpore.utils.glob.iglob")
def test_find_fastq_files(mock_iglob, mock_is_file):
    mock_is_file.return_value = True
    mock_iglob.return_value = [
        "./",
        "./CHANGELOG.md",
        "./tbpore",
        "./tbpore/constants.py",
        "./pyproject.toml",
        "./README.md",
        "./environment.yaml",
        "./tmp/toy/a.fq",
        "./tmp/toy/a.fq.gz",
        "./tmp/toy/a.fastq.gz",
        "./tmp/toy/a.fastq",
        "./tmp/toy/a.fasq.gz",
    ]

    actual = sorted(find_fastq_files(Path(".")))
    expected = sorted(
        [
            "./tmp/toy/a.fq",
            "./tmp/toy/a.fq.gz",
            "./tmp/toy/a.fastq.gz",
            "./tmp/toy/a.fastq",
        ]
    )

    assert actual == expected


class TestFastqPrefix:
    def test_no_fastq_prefix_returns_input_filename(self):
        path = "path/to/my.file"

        actual = fastq_prefix(path)
        expected = Path(path).name

        assert actual == expected

    def test_all_fastq_extensions_return_the_same(self):
        exts = ["fastq", "fq", "fastq.gz", "fq.gz"]

        expected = "my.file"
        for ext in exts:
            path = f"path/to/my.file.{ext}"

            actual = fastq_prefix(path)

            assert actual == expected


def test_download_file():
    url = "https://raw.githubusercontent.com/mbhall88/head_to_head_pipeline/c4e798608e9d9ffad5853ecda32007160feb500a/analysis/resistance_prediction/lsf.yaml"
    with tempfile.NamedTemporaryFile() as tmp:
        filename = Path(tmp.name)
        download_file(url, filename)
        result = filename.read_text()

    expected = """__default__:
  - "-E '${SOFTWAREDIR}/singularity_preexec_test/singularity_test.sh -t 20s /'"

"""

    assert result == expected


class TestValidateSha256:
    def test_matches(self):
        url = "https://raw.githubusercontent.com/mbhall88/head_to_head_pipeline/c4e798608e9d9ffad5853ecda32007160feb500a/analysis/resistance_prediction/lsf.yaml"
        expected_hash = (
            "50966ec6b51c8ad11d8e43e759b6fd0560b05590c90ceaed8b471ae448e68048"
        )
        with tempfile.NamedTemporaryFile() as tmp:
            filename = Path(tmp.name)
            download_file(url, filename)
            assert validate_sha256(filename, expected_hash)

    def test_doesnt_match(self):
        url = "https://raw.githubusercontent.com/mbhall88/head_to_head_pipeline/c4e798608e9d9ffad5853ecda32007160feb500a/analysis/resistance_prediction/lsf.yaml"
        expected_hash = (
            "40966ec6b51c8ad11d8e43e759b6fd0560b05590c90ceaed8b471ae448e68048"
        )
        with tempfile.NamedTemporaryFile() as tmp:
            filename = Path(tmp.name)
            download_file(url, filename)
            assert not validate_sha256(filename, expected_hash)


class TestDecompressFile:
    def test_decompress_and_dont_remove(self):
        p = Path("tmp.txt.gz")
        with gzip.open(p, "wb") as fp:
            fp.write("foo".encode())

        decompressed = p.with_suffix("")
        decompress_file(p, decompressed, remove_compressed=False)

        result = decompressed.read_text()
        expected = "foo"

        assert result == expected
        assert p.exists()
        p.unlink()
        decompressed.unlink()

    def test_decompress_and_remove(self):
        p = Path("tmp.txt.gz")
        with gzip.open(p, "wb") as fp:
            fp.write("foo".encode())

        decompressed = p.with_suffix("")
        decompress_file(p, decompressed, remove_compressed=True)

        result = decompressed.read_text()
        expected = "foo"

        assert result == expected
        assert not p.exists()
        decompressed.unlink()
