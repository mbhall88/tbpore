from pathlib import Path
from unittest import mock

from tbpore.utils import fastq_prefix, find_fastq_files, is_fastq


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
