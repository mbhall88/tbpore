"""Integration tests"""
import gzip
import sys

from click.testing import CliRunner
from unittest.mock import patch

from tbpore.constants import *
from tbpore.cli import main
from tbpore.external_tools import ExternalTool


@patch.object(ExternalTool, ExternalTool._run_core.__name__)
class TestCLIExecution:
    @staticmethod
    def get_command_line_from_mock(mock, index):
        return " ".join(mock.call_args_list[index].args[0])

    def test_whole_execution___minimum_params___check_all_external_tools_are_called_correctly(
            self, run_core_mock, tmp_path):
        runner = CliRunner()
        with runner.isolated_filesystem(temp_dir=tmp_path) as td:
            td = Path(td)
            infile = td / "in.fq"
            with open(infile, "w") as fp:
                fp.write("@r1\nACGT\n+$$$%\n")
            opts = [str(infile)]
            result = runner.invoke(main, opts)
            assert result.exit_code == 0

            # ensure all tools were called in the correct order and with the correct parameters
            assert run_core_mock.call_count == 8

            mykrobe_cl = self.get_command_line_from_mock(run_core_mock, 0)
            assert mykrobe_cl == \
                   f"mykrobe predict --sample in -t 1 --tmp .tbpore --skeleton_dir {cache_dir} -e 0.08 " \
                   f"--ploidy haploid --force --format json --min_proportion_expected_depth 0.20 --species tb " \
                   f"-m 8192MB -o .tbpore/in.mykrobe.json -i .tbpore/in.fq.gz"

            rasusa_cl = self.get_command_line_from_mock(run_core_mock, 1)
            assert rasusa_cl == \
                   f"rasusa -c 150 -g 4411532 -s 88 -o .tbpore/in.subsampled.fastq.gz -i .tbpore/in.fq.gz"

            minimap2_cl = self.get_command_line_from_mock(run_core_mock, 2)
            assert minimap2_cl == \
                   f"minimap2 -t 1 -a -L --sam-hit-only --secondary=no -x map-ont -o .tbpore/in.subsampled.sam " \
                   f"{H37RV_genome} .tbpore/in.subsampled.fastq.gz"

            samtools_sort_cl = self.get_command_line_from_mock(run_core_mock, 3)
            assert samtools_sort_cl == \
                   f"samtools sort -@ 1 -o .tbpore/in.subsampled.sorted.sam .tbpore/in.subsampled.sam"

            bcftools_mpileup_cl = self.get_command_line_from_mock(run_core_mock, 4)
            assert bcftools_mpileup_cl == \
                   f"bcftools mpileup -f {H37RV_genome} --threads 1 -x -O b -I -Q 13 " \
                   f"-a INFO/SCR,FORMAT/SP,INFO/ADR,INFO/ADF -h100 -M10000 -o .tbpore/in.subsampled.pileup.bcf " \
                   f".tbpore/in.subsampled.sorted.sam"

            bcftools_call_cl = self.get_command_line_from_mock(run_core_mock, 5)
            assert bcftools_call_cl == \
                   f"bcftools call --threads 1 --ploidy 1 -O b -V indels -m -o .tbpore/in.subsampled.snps.bcf " \
                   f".tbpore/in.subsampled.pileup.bcf"

            filter_vcf_cl = self.get_command_line_from_mock(run_core_mock, 6)
            assert filter_vcf_cl == \
                   f"{sys.executable} {external_scripts_dir}/apply_filters.py -P --verbose --overwrite -d 0 -D 0 " \
                   f"-q 85 -s 1 -b 0 -m 0 -r 0 -V 1e-05 -G 0 -K 0.9 -M 0 -x 0.2 " \
                   f"-o .tbpore/in.subsampled.snps.filtered.bcf -i .tbpore/in.subsampled.snps.bcf"

            generate_consensus_cl = self.get_command_line_from_mock(run_core_mock, 7)
            assert generate_consensus_cl == \
                   f"{sys.executable} {external_scripts_dir}/consensus.py --sample-id in --verbose --ignore all " \
                   f"--het-default none -o ./in.consensus.fa -i .tbpore/in.subsampled.snps.filtered.bcf " \
                   f"-f {H37RV_genome} -m {H37RV_mask}"


@patch.object(ExternalTool, ExternalTool._run_core.__name__)
class TestCLICleanup:
    def test_no_cleanup(self, run_core_mock, tmp_path):
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

    def test_with_cleanup(self, run_core_mock, tmp_path):
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


@patch.object(ExternalTool, ExternalTool._run_core.__name__)
class TestInputConcatenation:
    def test_single_file(self, run_core_mock, tmp_path):
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

    def test_multiple_files(self, run_core_mock, tmp_path):
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

    def test_input_is_dir(self, run_core_mock, tmp_path):
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

    def test_input_is_dir_and_one_file_has_bad_suffix(self, run_core_mock, tmp_path):
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
