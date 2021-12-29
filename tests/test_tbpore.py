"""Integration tests"""
import gzip
from pathlib import Path

from click.testing import CliRunner

from tbpore.tbpore import TMP_NAME, main


class TestCLICleanup:
    def test_no_cleanup(self, tmp_path):
        sample = "sam"
        opts = ["--no-cleanup", "-S", sample]
        runner = CliRunner()
        with runner.isolated_filesystem(temp_dir=tmp_path) as td:
            td = Path(td)
            infile = td / "in.fq"
            with open(infile, "w") as fp:
                fp.write("@r1\nACGT\n+$$$%\n")

            opts.extend(["-o", str(td)])
            opts.extend([str(infile)])
            result = runner.invoke(main, opts)
            assert result.exit_code == 0
            tbpore_tmp = td / TMP_NAME
            assert tbpore_tmp.exists()

    def test_with_cleanup(self, tmp_path):
        sample = "sam"
        opts = ["--cleanup", "-S", sample]
        runner = CliRunner()
        with runner.isolated_filesystem(temp_dir=tmp_path) as td:
            td = Path(td)
            infile = td / "in.fq"
            with open(infile, "w") as fp:
                fp.write("@r1\nACGT\n+$$$%\n")

            opts.extend(["-o", str(td)])
            opts.extend([str(infile)])
            result = runner.invoke(main, opts)
            assert result.exit_code == 0
            tbpore_tmp = td / TMP_NAME
            assert not tbpore_tmp.exists()


class TestInputConcatenation:
    def test_single_file(self, tmp_path):
        sample = "sam"
        opts = ["-D", "-S", sample]
        runner = CliRunner()
        expected_fq = "@r1\nACGT\n+$$$%\n"
        with runner.isolated_filesystem(temp_dir=tmp_path) as td:
            td = Path(td)
            infile = td / "in.fq"
            with open(infile, "w") as fp:
                fp.write(expected_fq)

            opts.extend(["-o", str(td)])
            opts.extend([str(infile)])
            result = runner.invoke(main, opts)
            assert result.exit_code == 0
            tbpore_concat = td / TMP_NAME / f"{sample}.fq.gz"
            assert tbpore_concat.exists()
            with gzip.open(tbpore_concat, mode="rt") as fp:
                actual = fp.read()
            assert actual == expected_fq

    def test_multiple_files(self, tmp_path):
        sample = "sam"
        opts = ["-D", "-S", sample]
        runner = CliRunner()
        expected_fq1 = "@r1\nACGT\n+$$$%\n"
        expected_fq2 = "@r1\nAAGT\n+$$$%\n"
        with runner.isolated_filesystem(temp_dir=tmp_path) as td:
            td = Path(td)
            infile1 = td / "in1.fq"
            infile2 = td / "in2.fastq.gz"
            with open(infile1, "w") as fp:
                fp.write(expected_fq1)
            with gzip.open(infile2, "wb") as fp:
                fp.write(expected_fq2.encode())

            opts.extend(["-o", str(td)])
            opts.extend([str(infile1), str(infile2)])
            result = runner.invoke(main, opts)
            assert result.exit_code == 0
            tbpore_concat = td / TMP_NAME / f"{sample}.fq.gz"
            assert tbpore_concat.exists()
            with gzip.open(tbpore_concat, mode="rt") as fp:
                actual = fp.read()

            expected = expected_fq1 + expected_fq2
            assert sorted(actual) == sorted(expected)

    def test_input_is_dir(self, tmp_path):
        sample = "sam"
        opts = ["-D", "-S", sample]
        runner = CliRunner()
        expected_fq1 = "@r1\nACGT\n+$$$%\n"
        expected_fq2 = "@r1\nAAGT\n+$$$%\n"
        with runner.isolated_filesystem(temp_dir=tmp_path) as td:
            td = Path(td)
            infile1 = td / "in1.fq"
            infile2 = td / "in2.fastq.gz"
            with open(infile1, "w") as fp:
                fp.write(expected_fq1)
            with gzip.open(infile2, "wb") as fp:
                fp.write(expected_fq2.encode())

            opts.extend(["-o", str(td)])
            opts.extend([str(td)])
            result = runner.invoke(main, opts)
            assert result.exit_code == 0
            tbpore_concat = td / TMP_NAME / f"{sample}.fq.gz"
            assert tbpore_concat.exists()
            with gzip.open(tbpore_concat, mode="rt") as fp:
                actual = fp.read()

            expected = expected_fq1 + expected_fq2
            assert sorted(actual) == sorted(expected)

    def test_input_is_dir_and_one_file_has_bad_suffix(self, tmp_path):
        sample = "sam"
        opts = ["-D", "-S", sample]
        runner = CliRunner()
        expected_fq1 = "@r1\nACGT\n+$$$%\n"
        expected_fq2 = "@r1\nAAGT\n+$$$%\n"
        with runner.isolated_filesystem(temp_dir=tmp_path) as td:
            td = Path(td)
            infile1 = td / "in1.faq"
            infile2 = td / "in2.fastq.gz"
            with open(infile1, "w") as fp:
                fp.write(expected_fq1)
            with gzip.open(infile2, "wb") as fp:
                fp.write(expected_fq2.encode())

            opts.extend(["-o", str(td)])
            opts.extend([str(td)])
            result = runner.invoke(main, opts)
            assert result.exit_code == 0
            tbpore_concat = td / TMP_NAME / f"{sample}.fq.gz"
            assert tbpore_concat.exists()
            with gzip.open(tbpore_concat, mode="rt") as fp:
                actual = fp.read()

            expected = expected_fq2
            assert sorted(actual) == sorted(expected)
