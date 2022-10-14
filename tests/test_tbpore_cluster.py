from pathlib import Path
from unittest.mock import patch

from click.testing import CliRunner

from tbpore.external_tools import ExternalTool
from tbpore.tbpore import TMP_NAME, main_cli


@patch.object(ExternalTool, ExternalTool._run_core.__name__)
class TestClusterCLI:
    @staticmethod
    def get_command_line_from_mock(mock, index):
        return " ".join(mock.call_args_list[index].args[0])

    def test_no_input___fails(self, run_core_mock, tmp_path):
        opts = ["cluster"]
        runner = CliRunner()
        with runner.isolated_filesystem(temp_dir=tmp_path):
            result = runner.invoke(main_cli, opts)
            assert result.exit_code == 2
            assert (
                b"To cluster consensus sequences, please provide at least two input consensus sequences"
                in result.stdout_bytes
            )

    def test_single_fasta_as_input___fails(self, run_core_mock, tmp_path):
        runner = CliRunner()
        with runner.isolated_filesystem(temp_dir=tmp_path) as td:
            td = Path(td)
            infile = td / "in.fq"
            with open(infile, "w") as fp:
                fp.write("@r1\nACGT\n+$$$%\n")
            opts = ["cluster", str(infile)]
            result = runner.invoke(main_cli, opts)
            assert result.exit_code == 2
            assert (
                b"To cluster consensus sequences, please provide at least two input consensus sequences"
                in result.stdout_bytes
            )

    @patch("tbpore.tbpore.produce_clusters")
    def test_whole_execution___minimum_params(
        self, produce_clusters_mock, run_core_mock, tmp_path
    ):
        runner = CliRunner()
        with runner.isolated_filesystem(temp_dir=tmp_path) as td:
            td = Path(td)
            infile_1 = td / "in1.fq"
            with open(infile_1, "w") as fp:
                fp.write("@r1\nACGT\n+$$$%\n")
            infile_2 = td / "in2.fq"
            with open(infile_2, "w") as fp:
                fp.write("@r1\nACGT\n+$$$%\n")
            opts = ["cluster", "-o", str(td), str(infile_1), str(infile_2)]
            result = runner.invoke(main_cli, opts)
            print(result)
            assert result.exit_code == 0

            assert run_core_mock.call_count == 1

            psdm_cl = self.get_command_line_from_mock(run_core_mock, 0)

            psdm_matrix = td / TMP_NAME / "psdm.matrix.csv"
            assert (
                psdm_cl
                == f"psdm --ignore-case --quiet --sort -t 1 -o {psdm_matrix} {td}/{TMP_NAME}/all_sequences.fq.gz"
            )

            threshold = 6
            produce_clusters_mock.assert_called_once_with(psdm_matrix, threshold, td)

    @patch("tbpore.tbpore.produce_clusters")
    def test_whole_execution___several_params(
        self, produce_clusters_mock, run_core_mock, tmp_path
    ):
        runner = CliRunner()
        with runner.isolated_filesystem(temp_dir=tmp_path) as td:
            td = Path(td)
            infile_1 = td / "in1.fq"
            with open(infile_1, "w") as fp:
                fp.write("@r1\nACGT\n+$$$%\n")
            infile_2 = td / "in2.fq"
            with open(infile_2, "w") as fp:
                fp.write("@r1\nACGT\n+$$$%\n")
            opts = [
                "cluster",
                "-o",
                str(td),
                "--threshold",
                "500",
                "--tmp",
                str(td / "custom_tmp"),
                "--threads",
                "101",
                "--cleanup",
                str(infile_1),
                str(infile_2),
            ]
            result = runner.invoke(main_cli, opts)
            assert result.exit_code == 0

            assert run_core_mock.call_count == 1

            psdm_cl = self.get_command_line_from_mock(run_core_mock, 0)

            psdm_matrix = td / "custom_tmp/psdm.matrix.csv"
            assert (
                psdm_cl
                == f"psdm --ignore-case --quiet --sort -t 101 -o {psdm_matrix} {td}/custom_tmp/all_sequences.fq.gz"
            )

            threshold = 500
            produce_clusters_mock.assert_called_once_with(psdm_matrix, threshold, td)

    @patch("tbpore.tbpore.produce_clusters")
    def test_no_cleanup(self, produce_clusters_mock, run_core_mock, tmp_path):
        runner = CliRunner()
        with runner.isolated_filesystem(temp_dir=tmp_path) as td:
            td = Path(td)
            infile_1 = td / "in1.fq"
            with open(infile_1, "w") as fp:
                fp.write("@r1\nACGT\n+$$$%\n")
            infile_2 = td / "in2.fq"
            with open(infile_2, "w") as fp:
                fp.write("@r1\nACGT\n+$$$%\n")
            opts = [
                "cluster",
                "--no-cleanup",
                "-o",
                str(td),
                str(infile_1),
                str(infile_2),
            ]
            result = runner.invoke(main_cli, opts)
            assert result.exit_code == 0

            tbpore_tmp = td / TMP_NAME
            assert tbpore_tmp.exists()

    @patch("tbpore.tbpore.produce_clusters")
    def test_with_cleanup(self, produce_clusters_mock, run_core_mock, tmp_path):
        runner = CliRunner()
        with runner.isolated_filesystem(temp_dir=tmp_path) as td:
            td = Path(td)
            infile_1 = td / "in1.fq"
            with open(infile_1, "w") as fp:
                fp.write("@r1\nACGT\n+$$$%\n")
            infile_2 = td / "in2.fq"
            with open(infile_2, "w") as fp:
                fp.write("@r1\nACGT\n+$$$%\n")
            opts = ["cluster", "--cleanup", "-o", str(td), str(infile_1), str(infile_2)]
            result = runner.invoke(main_cli, opts)
            assert result.exit_code == 0

            tbpore_tmp = td / TMP_NAME
            assert not tbpore_tmp.exists()
