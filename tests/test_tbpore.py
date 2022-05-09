"""Integration tests"""
import gzip
import subprocess
import sys
from pathlib import Path
from unittest.mock import patch

from click.testing import CliRunner

from tbpore.external_tools import ExternalTool
from tbpore.tbpore import (
    TMP_NAME,
    H37RV_genome,
    H37RV_mask,
    cache_dir,
    decontamination_db_fasta,
    decontamination_db_metadata,
    external_scripts_dir,
    main,
)


class TestExternalToolsExecution:
    @staticmethod
    def get_command_line_from_mock(mock, index):
        return " ".join(mock.call_args_list[index].args[0])

    @patch.object(ExternalTool, ExternalTool._run_core.__name__)
    def test_whole_execution___minimum_params___check_all_external_tools_are_called_correctly(
        self, run_core_mock, tmp_path
    ):
        runner = CliRunner()
        with runner.isolated_filesystem(temp_dir=tmp_path) as td:
            td = Path(td)
            infile = td / "in.fq"
            with open(infile, "w") as fp:
                fp.write("@r1\nACGT\n+$$$%\n")
            opts = [str(infile), "-o", str(td)]
            result = runner.invoke(main, opts)

            # check tbpore ran fine
            assert result.exit_code == 0

            # ensure all tools were called in the correct order and with the correct parameters
            assert run_core_mock.call_count == 14

            mykrobe_cl = self.get_command_line_from_mock(run_core_mock, 0)
            assert (
                mykrobe_cl
                == f"mykrobe predict --sample in -t 1 --tmp {td}/{TMP_NAME} --skeleton_dir {cache_dir} -e 0.08 "
                f"--ploidy haploid --format json --min_proportion_expected_depth 0.20 --species tb "
                f"-m 2048MB -o {td}/in.mykrobe.json -i {td}/{TMP_NAME}/in.fq.gz"
            )

            index_decontamination_db_cl = self.get_command_line_from_mock(
                run_core_mock, 1
            )
            assert (
                index_decontamination_db_cl
                == f"minimap2 -I 500M -x map-ont -t 1 -d {td}/{TMP_NAME}/tbpore.remove_contam.fa.gz.map-ont.mmi {decontamination_db_fasta}"
            )

            map_decontamination_db_cl = self.get_command_line_from_mock(
                run_core_mock, 2
            )
            assert (
                map_decontamination_db_cl
                == f"minimap2 -aL2 -x map-ont -t 1 -o {td}/{TMP_NAME}/in.decontaminated.sam {td}/{TMP_NAME}/tbpore.remove_contam.fa.gz.map-ont.mmi {td}/{TMP_NAME}/in.fq.gz"
            )

            sort_decontaminated_sam_cl = self.get_command_line_from_mock(
                run_core_mock, 3
            )
            assert (
                sort_decontaminated_sam_cl
                == f"samtools sort -@ 1 -o {td}/{TMP_NAME}/in.decontaminated.sorted.bam {td}/{TMP_NAME}/in.decontaminated.sam"
            )

            index_sorted_decontaminated_bam_cl = self.get_command_line_from_mock(
                run_core_mock, 4
            )
            assert (
                index_sorted_decontaminated_bam_cl
                == f"samtools index -@ 1 {td}/{TMP_NAME}/in.decontaminated.sorted.bam"
            )

            filter_contamination_cl = self.get_command_line_from_mock(run_core_mock, 5)
            assert (
                filter_contamination_cl
                == f"{sys.executable} {external_scripts_dir}/filter_contamination.py --verbose --ignore-secondary -o {td}/{TMP_NAME}/in.decontaminated.filter -i {td}/{TMP_NAME}/in.decontaminated.sorted.bam -m {decontamination_db_metadata}"
            )

            extract_decontaminated_nanopore_reads_cl = self.get_command_line_from_mock(
                run_core_mock, 6
            )
            assert (
                extract_decontaminated_nanopore_reads_cl
                == f"seqkit grep -o {td}/{TMP_NAME}/in.decontaminated.fastq.gz -f {td}/{TMP_NAME}/in.decontaminated.filter/keep.reads {td}/{TMP_NAME}/in.fq.gz"
            )

            rasusa_cl = self.get_command_line_from_mock(run_core_mock, 7)
            assert (
                rasusa_cl
                == f"rasusa -c 150 -g 4411532 -s 88 -o {td}/{TMP_NAME}/in.subsampled.fastq.gz -i {td}/{TMP_NAME}/in.decontaminated.fastq.gz"
            )

            minimap2_cl = self.get_command_line_from_mock(run_core_mock, 8)
            assert (
                minimap2_cl
                == f"minimap2 -t 1 -a -L --sam-hit-only --secondary=no -x map-ont -o {td}/{TMP_NAME}/in.subsampled.sam "
                f"{H37RV_genome} {td}/{TMP_NAME}/in.subsampled.fastq.gz"
            )

            samtools_sort_cl = self.get_command_line_from_mock(run_core_mock, 9)
            assert (
                samtools_sort_cl
                == f"samtools sort -@ 1 -o {td}/{TMP_NAME}/in.subsampled.sorted.sam {td}/{TMP_NAME}/in.subsampled.sam"
            )

            bcftools_mpileup_cl = self.get_command_line_from_mock(run_core_mock, 10)
            assert (
                bcftools_mpileup_cl
                == f"bcftools mpileup -f {H37RV_genome} --threads 1 -x -I -Q 13 "
                f"-a INFO/SCR,FORMAT/SP,INFO/ADR,INFO/ADF -h100 -M10000 -o {td}/{TMP_NAME}/in.subsampled.pileup.vcf "
                f"{td}/{TMP_NAME}/in.subsampled.sorted.sam"
            )

            bcftools_call_cl = self.get_command_line_from_mock(run_core_mock, 11)
            assert (
                bcftools_call_cl
                == f"bcftools call --threads 1 --ploidy 1 -V indels -m -o {td}/{TMP_NAME}/in.subsampled.snps.vcf "
                f"{td}/{TMP_NAME}/in.subsampled.pileup.vcf"
            )

            filter_vcf_cl = self.get_command_line_from_mock(run_core_mock, 12)
            assert (
                filter_vcf_cl
                == f"{sys.executable} {external_scripts_dir}/apply_filters.py -P --verbose --overwrite -d 0 -D 0 "
                f"-q 85 -s 1 -b 0 -m 0 -r 0 -V 1e-05 -G 0 -K 0.9 -M 0 -x 0.2 "
                f"-o {td}/in.subsampled.snps.filtered.bcf -i {td}/{TMP_NAME}/in.subsampled.snps.vcf"
            )

            generate_consensus_cl = self.get_command_line_from_mock(run_core_mock, 13)
            assert (
                generate_consensus_cl
                == f"{sys.executable} {external_scripts_dir}/consensus.py --sample-id in --verbose --ignore all "
                f"--het-default none -o {td}/in.consensus.fa -i {td}/in.subsampled.snps.filtered.bcf "
                f"-f {H37RV_genome} -m {H37RV_mask}"
            )

    @patch.object(ExternalTool, ExternalTool._run_core.__name__)
    def test_whole_execution___several_params_affecting_tools_command_lines___check_all_external_tools_are_called_correctly(
        self, run_core_mock, tmp_path
    ):
        runner = CliRunner()
        with runner.isolated_filesystem(temp_dir=tmp_path) as td:
            td = Path(td)
            infile = td / "in.fq"
            with open(infile, "w") as fp:
                fp.write("@r1\nACGT\n+$$$%\n")
            opts = [
                str(infile),
                "-o",
                str(td),
                "--tmp",
                str(td / "custom_tmp"),
                "--name",
                "custom_name",
                "--threads",
                "8",
                "--report_all_mykrobe_calls",
            ]
            result = runner.invoke(main, opts)

            # check tbpore ran fine
            assert result.exit_code == 0

            # ensure all tools were called in the correct order and with the correct parameters
            assert run_core_mock.call_count == 14

            mykrobe_cl = self.get_command_line_from_mock(run_core_mock, 0)
            assert (
                mykrobe_cl
                == f"mykrobe predict -A --sample custom_name -t 8 --tmp {td}/custom_tmp --skeleton_dir {cache_dir} -e 0.08 "
                f"--ploidy haploid --format json --min_proportion_expected_depth 0.20 --species tb "
                f"-m 2048MB -o {td}/custom_name.mykrobe.json -i {td}/custom_tmp/custom_name.fq.gz"
            )

            index_decontamination_db_cl = self.get_command_line_from_mock(
                run_core_mock, 1
            )
            assert (
                index_decontamination_db_cl
                == f"minimap2 -I 500M -x map-ont -t 8 -d {td}/custom_tmp/tbpore.remove_contam.fa.gz.map-ont.mmi {decontamination_db_fasta}"
            )

            map_decontamination_db_cl = self.get_command_line_from_mock(
                run_core_mock, 2
            )
            assert (
                map_decontamination_db_cl
                == f"minimap2 -aL2 -x map-ont -t 8 -o {td}/custom_tmp/custom_name.decontaminated.sam {td}/custom_tmp/tbpore.remove_contam.fa.gz.map-ont.mmi {td}/custom_tmp/custom_name.fq.gz"
            )

            sort_decontaminated_sam_cl = self.get_command_line_from_mock(
                run_core_mock, 3
            )
            assert (
                sort_decontaminated_sam_cl
                == f"samtools sort -@ 8 -o {td}/custom_tmp/custom_name.decontaminated.sorted.bam {td}/custom_tmp/custom_name.decontaminated.sam"
            )

            index_sorted_decontaminated_bam_cl = self.get_command_line_from_mock(
                run_core_mock, 4
            )
            assert (
                index_sorted_decontaminated_bam_cl
                == f"samtools index -@ 8 {td}/custom_tmp/custom_name.decontaminated.sorted.bam"
            )

            filter_contamination_cl = self.get_command_line_from_mock(run_core_mock, 5)
            assert (
                filter_contamination_cl
                == f"{sys.executable} {external_scripts_dir}/filter_contamination.py --verbose --ignore-secondary -o {td}/custom_tmp/custom_name.decontaminated.filter -i {td}/custom_tmp/custom_name.decontaminated.sorted.bam -m {decontamination_db_metadata}"
            )

            extract_decontaminated_nanopore_reads_cl = self.get_command_line_from_mock(
                run_core_mock, 6
            )
            assert (
                extract_decontaminated_nanopore_reads_cl
                == f"seqkit grep -o {td}/custom_tmp/custom_name.decontaminated.fastq.gz -f {td}/custom_tmp/custom_name.decontaminated.filter/keep.reads {td}/custom_tmp/custom_name.fq.gz"
            )

            rasusa_cl = self.get_command_line_from_mock(run_core_mock, 7)
            assert (
                rasusa_cl
                == f"rasusa -c 150 -g 4411532 -s 88 -o {td}/custom_tmp/custom_name.subsampled.fastq.gz "
                f"-i {td}/custom_tmp/custom_name.decontaminated.fastq.gz"
            )

            minimap2_cl = self.get_command_line_from_mock(run_core_mock, 8)
            assert (
                minimap2_cl
                == f"minimap2 -t 8 -a -L --sam-hit-only --secondary=no -x map-ont -o {td}/custom_tmp/custom_name.subsampled.sam "
                f"{H37RV_genome} {td}/custom_tmp/custom_name.subsampled.fastq.gz"
            )

            samtools_sort_cl = self.get_command_line_from_mock(run_core_mock, 9)
            assert (
                samtools_sort_cl
                == f"samtools sort -@ 8 -o {td}/custom_tmp/custom_name.subsampled.sorted.sam {td}/custom_tmp/custom_name.subsampled.sam"
            )

            bcftools_mpileup_cl = self.get_command_line_from_mock(run_core_mock, 10)
            assert (
                bcftools_mpileup_cl
                == f"bcftools mpileup -f {H37RV_genome} --threads 8 -x -I -Q 13 "
                f"-a INFO/SCR,FORMAT/SP,INFO/ADR,INFO/ADF -h100 -M10000 -o {td}/custom_tmp/custom_name.subsampled.pileup.vcf "
                f"{td}/custom_tmp/custom_name.subsampled.sorted.sam"
            )

            bcftools_call_cl = self.get_command_line_from_mock(run_core_mock, 11)
            assert (
                bcftools_call_cl
                == f"bcftools call --threads 8 --ploidy 1 -V indels -m -o {td}/custom_tmp/custom_name.subsampled.snps.vcf "
                f"{td}/custom_tmp/custom_name.subsampled.pileup.vcf"
            )

            filter_vcf_cl = self.get_command_line_from_mock(run_core_mock, 12)
            assert (
                filter_vcf_cl
                == f"{sys.executable} {external_scripts_dir}/apply_filters.py -P --verbose --overwrite -d 0 -D 0 "
                f"-q 85 -s 1 -b 0 -m 0 -r 0 -V 1e-05 -G 0 -K 0.9 -M 0 -x 0.2 "
                f"-o {td}/custom_name.subsampled.snps.filtered.bcf -i {td}/custom_tmp/custom_name.subsampled.snps.vcf"
            )

            generate_consensus_cl = self.get_command_line_from_mock(run_core_mock, 13)
            assert (
                generate_consensus_cl
                == f"{sys.executable} {external_scripts_dir}/consensus.py --sample-id custom_name --verbose --ignore all "
                f"--het-default none -o {td}/custom_name.consensus.fa -i {td}/custom_name.subsampled.snps.filtered.bcf "
                f"-f {H37RV_genome} -m {H37RV_mask}"
            )

    @patch.object(
        ExternalTool,
        ExternalTool._run_core.__name__,
        side_effect=["", "", subprocess.CalledProcessError(1, "minimap2")],
    )
    def test_partial_execution___minimum_params___minimap2_map_decontamination_db_fails___checks_fail_happens_and_previous_tools_called_correctly(
        self, run_core_mock, tmp_path
    ):
        runner = CliRunner()
        with runner.isolated_filesystem(temp_dir=tmp_path) as td:
            td = Path(td)
            infile = td / "in.fq"
            with open(infile, "w") as fp:
                fp.write("@r1\nACGT\n+$$$%\n")
            opts = [
                str(infile),
                "--cleanup",
                "-o",
                str(td),
            ]  # asks for cleanup, but it should not cleanup as this is a run that fails
            result = runner.invoke(main, opts)

            # check if tbpore indeed failed
            assert result.exit_code == 1
            map_decontamination_db_cl = self.get_command_line_from_mock(
                run_core_mock, 2
            )
            assert (
                b"Error calling "
                + map_decontamination_db_cl.encode("utf-8")
                + b" (return code 1)"
                in result.stdout_bytes
            )

            # check if all tools until minimap2 were called correctly
            assert run_core_mock.call_count == 3

            mykrobe_cl = self.get_command_line_from_mock(run_core_mock, 0)
            assert (
                mykrobe_cl
                == f"mykrobe predict --sample in -t 1 --tmp {td}/{TMP_NAME} --skeleton_dir {cache_dir} -e 0.08 "
                f"--ploidy haploid --format json --min_proportion_expected_depth 0.20 --species tb "
                f"-m 2048MB -o {td}/in.mykrobe.json -i {td}/{TMP_NAME}/in.fq.gz"
            )

            index_decontamination_db_cl = self.get_command_line_from_mock(
                run_core_mock, 1
            )
            assert (
                index_decontamination_db_cl
                == f"minimap2 -I 500M -x map-ont -t 1 -d {td}/{TMP_NAME}/tbpore.remove_contam.fa.gz.map-ont.mmi {decontamination_db_fasta}"
            )

            assert (
                map_decontamination_db_cl
                == f"minimap2 -aL2 -x map-ont -t 1 -o {td}/{TMP_NAME}/in.decontaminated.sam {td}/{TMP_NAME}/tbpore.remove_contam.fa.gz.map-ont.mmi {td}/{TMP_NAME}/in.fq.gz"
            )

            # check if tmp not removed
            tbpore_tmp = td / TMP_NAME
            assert tbpore_tmp.exists()


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
    def test_no_input___fails(self, run_core_mock, tmp_path):
        runner = CliRunner()
        with runner.isolated_filesystem(temp_dir=tmp_path):
            result = runner.invoke(main, [])
            assert result.exit_code == 2
            assert b"No INPUT files given" in result.stdout_bytes

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

    def test_single_gz_file(self, run_core_mock, tmp_path):
        sample = "sam"
        opts = ["-D", "-S", sample]
        runner = CliRunner()
        expected_fq = "@r1\nACGT\n+$$$%\n"
        with runner.isolated_filesystem(temp_dir=tmp_path) as td:
            td = Path(td)
            infile = td / "in.fastq.gz"
            with gzip.open(infile, "wb") as fp:
                fp.write(expected_fq.encode())

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

    def test_input_is_empty_dir___error_out(self, run_core_mock, tmp_path):
        sample = "sam"
        opts = ["-D", "-S", sample]
        runner = CliRunner()
        with runner.isolated_filesystem(temp_dir=tmp_path) as td:
            td = Path(td)
            opts.extend(["-o", str(td)])
            opts.extend([str(td)])
            result = runner.invoke(main, opts)
            assert result.exit_code == 2
            assert (
                b"No fastq files found for the given inputs, please check your input files/dirs."
                in result.stdout_bytes
            )

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

    def test_input_is_dir_and_files_in_dir___ensure_duplication_does_not_happen(
        self, run_core_mock, tmp_path
    ):
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
            opts.extend([str(td), str(infile1), str(infile2)])
            result = runner.invoke(main, opts)
            assert result.exit_code == 0
            print(result.stdout_bytes)
            assert b"Found 2 fastq files. Joining them..." in result.stdout_bytes

            tbpore_concat = td / TMP_NAME / f"{sample}.fq.gz"
            assert tbpore_concat.exists()
            with gzip.open(tbpore_concat, mode="rt") as fp:
                actual = fp.read()

            expected = expected_fq1 + expected_fq2
            assert sorted(actual) == sorted(expected)
