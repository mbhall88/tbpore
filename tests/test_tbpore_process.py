"""Integration tests"""
import gzip
import subprocess
import sys
from pathlib import Path
from unittest.mock import patch

from click.testing import CliRunner

from tbpore.external_tools import ExternalTool
from tbpore.tbpore import (
    CACHE_DIR,
    DECONTAMINATION_DB_INDEX,
    DECONTAMINATION_DB_METADATA,
    EXTERNAL_SCRIPTS_DIR,
    TMP_NAME,
    H37RV_genome,
    H37RV_mask,
    main_cli,
)


@patch("tbpore.tbpore.ensure_decontamination_db_is_available")
class TestExternalToolsExecution:
    @staticmethod
    def get_command_line_from_mock(mock, index):
        return " ".join(mock.call_args_list[index].args[0])

    @patch.object(ExternalTool, ExternalTool._run_core.__name__)
    def test_whole_execution___minimum_params___check_all_external_tools_are_called_correctly(
        self, run_core_mock, ensure_decontamination_db_is_available_mock, tmp_path
    ):
        runner = CliRunner()
        with runner.isolated_filesystem(temp_dir=tmp_path) as td:
            td = Path(td)
            infile = td / "in.fq"
            with open(infile, "w") as fp:
                fp.write("@r1\nACGT\n+$$$%\n")
            opts = ["process", str(infile), "-o", str(td)]
            result = runner.invoke(main_cli, opts)

            # check tbpore ran fine
            assert result.exit_code == 0, result.stderr

            # ensure all tools were called in the correct order and with the correct parameters
            assert run_core_mock.call_count == 15

            map_decontamination_db_cl = self.get_command_line_from_mock(
                run_core_mock, 0
            )
            assert (
                map_decontamination_db_cl
                == f"minimap2 -aL2 -x map-ont -t 1 -o {td}/{TMP_NAME}/in.decontaminated.sam {DECONTAMINATION_DB_INDEX} {td}/{TMP_NAME}/in.fq.gz"
            )

            sort_decontaminated_sam_cl = self.get_command_line_from_mock(
                run_core_mock, 1
            )
            assert (
                sort_decontaminated_sam_cl
                == f"samtools sort -@ 1 -o {td}/{TMP_NAME}/in.decontaminated.sorted.bam {td}/{TMP_NAME}/in.decontaminated.sam"
            )

            index_sorted_decontaminated_bam_cl = self.get_command_line_from_mock(
                run_core_mock, 2
            )
            assert (
                index_sorted_decontaminated_bam_cl
                == f"samtools index -@ 1 {td}/{TMP_NAME}/in.decontaminated.sorted.bam"
            )

            filter_contamination_cl = self.get_command_line_from_mock(run_core_mock, 3)
            assert (
                filter_contamination_cl
                == f"{sys.executable} {EXTERNAL_SCRIPTS_DIR}/filter_contamination.py --verbose --ignore-secondary -o {td}/{TMP_NAME}/in.decontaminated.filter -i {td}/{TMP_NAME}/in.decontaminated.sorted.bam -m {DECONTAMINATION_DB_METADATA}"
            )

            extract_decontaminated_nanopore_reads_cl = self.get_command_line_from_mock(
                run_core_mock, 4
            )
            assert (
                extract_decontaminated_nanopore_reads_cl
                == f"seqkit grep -o {td}/{TMP_NAME}/in.decontaminated.fastq.gz -f {td}/{TMP_NAME}/in.decontaminated.filter/keep.reads {td}/{TMP_NAME}/in.fq.gz"
            )

            sort_decontaminated_reads_cl = self.get_command_line_from_mock(
                run_core_mock, 5
            )
            assert (
                sort_decontaminated_reads_cl
                == f"seqkit sort -o {td}/{TMP_NAME}/in.sorted.fastq.gz {td}/{TMP_NAME}/in.decontaminated.fastq.gz"
            )

            rasusa_cl = self.get_command_line_from_mock(run_core_mock, 6)
            assert (
                rasusa_cl
                == f"rasusa -c 150 -g 4411532 -s 88 -o {td}/{TMP_NAME}/in.subsampled.fastq.gz -i {td}/{TMP_NAME}/in.sorted.fastq.gz"
            )

            nanoq_cl = self.get_command_line_from_mock(run_core_mock, 7)
            assert (
                nanoq_cl
                == f"nanoq -vv -s -r {td}/in.stats.txt -i {td}/{TMP_NAME}/in.subsampled.fastq.gz"
            )

            mykrobe_cl = self.get_command_line_from_mock(run_core_mock, 8)
            assert (
                mykrobe_cl
                == f"mykrobe predict --sample in -t 1 --tmp {td}/{TMP_NAME} --skeleton_dir {CACHE_DIR} --ont "
                f"--format json --min_proportion_expected_depth 0.20 --species tb "
                f"-m 2048MB -o {td}/in.mykrobe.json -i {td}/{TMP_NAME}/in.subsampled.fastq.gz"
            )

            minimap2_cl = self.get_command_line_from_mock(run_core_mock, 9)
            assert (
                minimap2_cl
                == f"minimap2 -t 1 -a -L --sam-hit-only --secondary=no -x map-ont -o {td}/{TMP_NAME}/in.subsampled.sam "
                f"{H37RV_genome} {td}/{TMP_NAME}/in.subsampled.fastq.gz"
            )

            samtools_sort_cl = self.get_command_line_from_mock(run_core_mock, 10)
            assert (
                samtools_sort_cl
                == f"samtools sort -@ 1 -o {td}/{TMP_NAME}/in.subsampled.sorted.sam {td}/{TMP_NAME}/in.subsampled.sam"
            )

            bcftools_mpileup_cl = self.get_command_line_from_mock(run_core_mock, 11)
            assert (
                bcftools_mpileup_cl
                == f"bcftools mpileup -f {H37RV_genome} --threads 1 -x -I -Q 13 "
                f"-a INFO/SCR,FORMAT/SP,INFO/ADR,INFO/ADF -h100 -M10000 -o {td}/{TMP_NAME}/in.pileup.vcf "
                f"{td}/{TMP_NAME}/in.subsampled.sorted.sam"
            )

            bcftools_call_cl = self.get_command_line_from_mock(run_core_mock, 12)
            assert (
                bcftools_call_cl
                == f"bcftools call --threads 1 --ploidy 1 -V indels -m -o {td}/{TMP_NAME}/in.snps.vcf "
                f"{td}/{TMP_NAME}/in.pileup.vcf"
            )

            filter_vcf_cl = self.get_command_line_from_mock(run_core_mock, 13)
            assert (
                filter_vcf_cl
                == f"{sys.executable} {EXTERNAL_SCRIPTS_DIR}/apply_filters.py -P --verbose --overwrite -d 5 -D 0 "
                f"-q 25 -s 1 -b 0 -m 0 -r 0 -V 1e-05 -G 0 -K 0.9 -M 30 -x 0.2 "
                f"-o {td}/in.snps.filtered.bcf -i {td}/{TMP_NAME}/in.snps.vcf"
            )

            generate_consensus_cl = self.get_command_line_from_mock(run_core_mock, 14)
            assert (
                generate_consensus_cl
                == f"{sys.executable} {EXTERNAL_SCRIPTS_DIR}/consensus.py --sample-id in --verbose --ignore all "
                f"--het-default none -o {td}/in.consensus.fa -i {td}/in.snps.filtered.bcf "
                f"-f {H37RV_genome} -m {H37RV_mask}"
            )

    @patch.object(ExternalTool, ExternalTool._run_core.__name__)
    def test_whole_execution___minimum_params___check_coverage_zero_disables_subsampling(
        self, run_core_mock, ensure_decontamination_db_is_available_mock, tmp_path
    ):
        runner = CliRunner()
        with runner.isolated_filesystem(temp_dir=tmp_path) as td:
            td = Path(td)
            infile = td / "in.fq"
            with open(infile, "w") as fp:
                fp.write("@r1\nACGT\n+$$$%\n")
            opts = ["process", str(infile), "-o", str(td), "--coverage", "0"]
            result = runner.invoke(main_cli, opts)

            # check tbpore ran fine
            assert result.exit_code == 0, result.stderr

            # ensure all tools were called in the correct order and with the correct parameters
            assert run_core_mock.call_count == 14

            map_decontamination_db_cl = self.get_command_line_from_mock(
                run_core_mock, 0
            )
            assert (
                map_decontamination_db_cl
                == f"minimap2 -aL2 -x map-ont -t 1 -o {td}/{TMP_NAME}/in.decontaminated.sam {DECONTAMINATION_DB_INDEX} {td}/{TMP_NAME}/in.fq.gz"
            )

            sort_decontaminated_sam_cl = self.get_command_line_from_mock(
                run_core_mock, 1
            )
            assert (
                sort_decontaminated_sam_cl
                == f"samtools sort -@ 1 -o {td}/{TMP_NAME}/in.decontaminated.sorted.bam {td}/{TMP_NAME}/in.decontaminated.sam"
            )

            index_sorted_decontaminated_bam_cl = self.get_command_line_from_mock(
                run_core_mock, 2
            )
            assert (
                index_sorted_decontaminated_bam_cl
                == f"samtools index -@ 1 {td}/{TMP_NAME}/in.decontaminated.sorted.bam"
            )

            filter_contamination_cl = self.get_command_line_from_mock(run_core_mock, 3)
            assert (
                filter_contamination_cl
                == f"{sys.executable} {EXTERNAL_SCRIPTS_DIR}/filter_contamination.py --verbose --ignore-secondary -o {td}/{TMP_NAME}/in.decontaminated.filter -i {td}/{TMP_NAME}/in.decontaminated.sorted.bam -m {DECONTAMINATION_DB_METADATA}"
            )

            extract_decontaminated_nanopore_reads_cl = self.get_command_line_from_mock(
                run_core_mock, 4
            )
            assert (
                extract_decontaminated_nanopore_reads_cl
                == f"seqkit grep -o {td}/{TMP_NAME}/in.decontaminated.fastq.gz -f {td}/{TMP_NAME}/in.decontaminated.filter/keep.reads {td}/{TMP_NAME}/in.fq.gz"
            )

            sort_decontaminated_reads_cl = self.get_command_line_from_mock(
                run_core_mock, 5
            )
            assert (
                sort_decontaminated_reads_cl
                == f"seqkit sort -o {td}/{TMP_NAME}/in.sorted.fastq.gz {td}/{TMP_NAME}/in.decontaminated.fastq.gz"
            )

            nanoq_cl = self.get_command_line_from_mock(run_core_mock, 6)
            assert (
                nanoq_cl
                == f"nanoq -vv -s -r {td}/in.stats.txt -i {td}/{TMP_NAME}/in.sorted.fastq.gz"
            )

            mykrobe_cl = self.get_command_line_from_mock(run_core_mock, 7)
            assert (
                mykrobe_cl
                == f"mykrobe predict --sample in -t 1 --tmp {td}/{TMP_NAME} --skeleton_dir {CACHE_DIR} --ont "
                f"--format json --min_proportion_expected_depth 0.20 --species tb "
                f"-m 2048MB -o {td}/in.mykrobe.json -i {td}/{TMP_NAME}/in.sorted.fastq.gz"
            )

            minimap2_cl = self.get_command_line_from_mock(run_core_mock, 8)
            assert (
                minimap2_cl
                == f"minimap2 -t 1 -a -L --sam-hit-only --secondary=no -x map-ont -o {td}/{TMP_NAME}/in.sam "
                f"{H37RV_genome} {td}/{TMP_NAME}/in.sorted.fastq.gz"
            )

            samtools_sort_cl = self.get_command_line_from_mock(run_core_mock, 9)
            assert (
                samtools_sort_cl
                == f"samtools sort -@ 1 -o {td}/{TMP_NAME}/in.sorted.sam {td}/{TMP_NAME}/in.sam"
            )

            bcftools_mpileup_cl = self.get_command_line_from_mock(run_core_mock, 10)
            assert (
                bcftools_mpileup_cl
                == f"bcftools mpileup -f {H37RV_genome} --threads 1 -x -I -Q 13 "
                f"-a INFO/SCR,FORMAT/SP,INFO/ADR,INFO/ADF -h100 -M10000 -o {td}/{TMP_NAME}/in.pileup.vcf "
                f"{td}/{TMP_NAME}/in.sorted.sam"
            )

            bcftools_call_cl = self.get_command_line_from_mock(run_core_mock, 11)
            assert (
                bcftools_call_cl
                == f"bcftools call --threads 1 --ploidy 1 -V indels -m -o {td}/{TMP_NAME}/in.snps.vcf "
                f"{td}/{TMP_NAME}/in.pileup.vcf"
            )

            filter_vcf_cl = self.get_command_line_from_mock(run_core_mock, 12)
            assert (
                filter_vcf_cl
                == f"{sys.executable} {EXTERNAL_SCRIPTS_DIR}/apply_filters.py -P --verbose --overwrite -d 5 -D 0 "
                f"-q 25 -s 1 -b 0 -m 0 -r 0 -V 1e-05 -G 0 -K 0.9 -M 30 -x 0.2 "
                f"-o {td}/in.snps.filtered.bcf -i {td}/{TMP_NAME}/in.snps.vcf"
            )

            generate_consensus_cl = self.get_command_line_from_mock(run_core_mock, 13)
            assert (
                generate_consensus_cl
                == f"{sys.executable} {EXTERNAL_SCRIPTS_DIR}/consensus.py --sample-id in --verbose --ignore all "
                f"--het-default none -o {td}/in.consensus.fa -i {td}/in.snps.filtered.bcf "
                f"-f {H37RV_genome} -m {H37RV_mask}"
            )

    @patch.object(ExternalTool, ExternalTool._run_core.__name__)
    def test_whole_execution___minimum_params___check_mapping_stats_added(
        self, run_core_mock, ensure_decontamination_db_is_available_mock, tmp_path
    ):
        runner = CliRunner()
        with runner.isolated_filesystem(temp_dir=tmp_path) as td:
            td = Path(td)
            infile = td / "in.fq"
            with open(infile, "w") as fp:
                fp.write("@r1\nACGT\n+$$$%\n")

            stats_report = Path(f"{td}/in.stats.txt")
            contam_dir = Path(f"{td}/{TMP_NAME}/in.decontaminated.filter")
            contam_dir.mkdir(parents=True, exist_ok=True)
            n_keep = 10000
            with open(contam_dir / "keep.reads", "w") as fp:
                for _ in range(n_keep):
                    print("foo", file=fp)

            n_contam = 412
            with open(contam_dir / "contaminant.reads", "w") as fp:
                for _ in range(n_contam):
                    print("foo", file=fp)

            n_unmapped = 22
            with open(contam_dir / "unmapped.reads", "w") as fp:
                for _ in range(n_unmapped):
                    print("foo", file=fp)

            original_stats = """Nanoq Read Summary
====================

Number of reads:      158583
Number of bases:      191101235
N50 read length:      1615
Longest read:         45649
Shortest read:        76
Mean read length:     1205
Median read length:   833
Mean read quality:    11.66
Median read quality:  11.95


Read length thresholds (bp)

> 200       158158            99.7%
> 500       122965            77.5%
> 1000      64833             40.9%
> 2000      22688             14.3%
> 5000      2603              01.6%
> 10000     325               00.2%
> 30000     2                 00.0%
> 50000     0                 00.0%
> 100000    0                 00.0%
> 1000000   0                 00.0%


Read quality thresholds (Q)

> 5   158506        100.0%
> 7   154884        97.7%
> 10  116393        73.4%
> 12  77982         49.2%
> 15  10451         06.6%
> 20  20            00.0%
> 25  0             00.0%
> 30  0             00.0%
"""
            stats_report.write_text(original_stats)

            opts = ["process", str(infile), "-o", str(td)]
            result = runner.invoke(main_cli, opts)

            # check tbpore ran fine
            assert result.exit_code == 0

            expected_stats = """Nanoq Read Summary
====================

Number of reads:      158583
Num. MTB reads:       10000 (95.84%)
Num. contam. reads:   412 (3.95%)
Num. unmapped reads:  22 (0.21%)
Number of bases:      191101235
N50 read length:      1615
Longest read:         45649
Shortest read:        76
Mean read length:     1205
Median read length:   833
Mean read quality:    11.66
Median read quality:  11.95


Read length thresholds (bp)

> 200       158158            99.7%
> 500       122965            77.5%
> 1000      64833             40.9%
> 2000      22688             14.3%
> 5000      2603              01.6%
> 10000     325               00.2%
> 30000     2                 00.0%
> 50000     0                 00.0%
> 100000    0                 00.0%
> 1000000   0                 00.0%


Read quality thresholds (Q)

> 5   158506        100.0%
> 7   154884        97.7%
> 10  116393        73.4%
> 12  77982         49.2%
> 15  10451         06.6%
> 20  20            00.0%
> 25  0             00.0%
> 30  0             00.0%
""".strip()
            actual_stats = stats_report.read_text().strip()

            assert actual_stats == expected_stats, print(actual_stats)

    @patch.object(ExternalTool, ExternalTool._run_core.__name__)
    def test_whole_execution___several_params_affecting_tools_command_lines___check_all_external_tools_are_called_correctly(
        self, run_core_mock, ensure_decontamination_db_is_available_mock, tmp_path
    ):
        runner = CliRunner()
        with runner.isolated_filesystem(temp_dir=tmp_path) as td:
            td = Path(td)
            infile = td / "in.fq"
            with open(infile, "w") as fp:
                fp.write("@r1\nACGT\n+$$$%\n")
            opts = [
                "process",
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
            result = runner.invoke(main_cli, opts)

            # check tbpore ran fine
            assert result.exit_code == 0, result.stderr

            # ensure all tools were called in the correct order and with the correct parameters
            assert run_core_mock.call_count == 15

            map_decontamination_db_cl = self.get_command_line_from_mock(
                run_core_mock, 0
            )
            assert (
                map_decontamination_db_cl
                == f"minimap2 -aL2 -x map-ont -t 8 -o {td}/custom_tmp/custom_name.decontaminated.sam {DECONTAMINATION_DB_INDEX} {td}/custom_tmp/custom_name.fq.gz"
            )

            sort_decontaminated_sam_cl = self.get_command_line_from_mock(
                run_core_mock, 1
            )
            assert (
                sort_decontaminated_sam_cl
                == f"samtools sort -@ 8 -o {td}/custom_tmp/custom_name.decontaminated.sorted.bam {td}/custom_tmp/custom_name.decontaminated.sam"
            )

            index_sorted_decontaminated_bam_cl = self.get_command_line_from_mock(
                run_core_mock, 2
            )
            assert (
                index_sorted_decontaminated_bam_cl
                == f"samtools index -@ 8 {td}/custom_tmp/custom_name.decontaminated.sorted.bam"
            )

            filter_contamination_cl = self.get_command_line_from_mock(run_core_mock, 3)
            assert (
                filter_contamination_cl
                == f"{sys.executable} {EXTERNAL_SCRIPTS_DIR}/filter_contamination.py --verbose --ignore-secondary -o {td}/custom_tmp/custom_name.decontaminated.filter -i {td}/custom_tmp/custom_name.decontaminated.sorted.bam -m {DECONTAMINATION_DB_METADATA}"
            )

            extract_decontaminated_nanopore_reads_cl = self.get_command_line_from_mock(
                run_core_mock, 4
            )
            assert (
                extract_decontaminated_nanopore_reads_cl
                == f"seqkit grep -o {td}/custom_tmp/custom_name.decontaminated.fastq.gz -f {td}/custom_tmp/custom_name.decontaminated.filter/keep.reads {td}/custom_tmp/custom_name.fq.gz"
            )

            sort_decontaminated_reads_cl = self.get_command_line_from_mock(
                run_core_mock, 5
            )
            assert (
                sort_decontaminated_reads_cl
                == f"seqkit sort -o {td}/custom_tmp/custom_name.sorted.fastq.gz {td}/custom_tmp/custom_name.decontaminated.fastq.gz"
            )

            rasusa_cl = self.get_command_line_from_mock(run_core_mock, 6)
            assert (
                rasusa_cl
                == f"rasusa -c 150 -g 4411532 -s 88 -o {td}/custom_tmp/custom_name.subsampled.fastq.gz "
                f"-i {td}/custom_tmp/custom_name.sorted.fastq.gz"
            )

            nanoq_cl = self.get_command_line_from_mock(run_core_mock, 7)
            assert (
                nanoq_cl
                == f"nanoq -vv -s -r {td}/custom_name.stats.txt -i {td}/custom_tmp/custom_name.subsampled.fastq.gz"
            )

            mykrobe_cl = self.get_command_line_from_mock(run_core_mock, 8)
            assert (
                mykrobe_cl
                == f"mykrobe predict -A --sample custom_name -t 8 --tmp {td}/custom_tmp --skeleton_dir {CACHE_DIR} --ont "
                f"--format json --min_proportion_expected_depth 0.20 --species tb "
                f"-m 2048MB -o {td}/custom_name.mykrobe.json -i {td}/custom_tmp/custom_name.subsampled.fastq.gz"
            )

            minimap2_cl = self.get_command_line_from_mock(run_core_mock, 9)
            assert (
                minimap2_cl
                == f"minimap2 -t 8 -a -L --sam-hit-only --secondary=no -x map-ont -o {td}/custom_tmp/custom_name.subsampled.sam "
                f"{H37RV_genome} {td}/custom_tmp/custom_name.subsampled.fastq.gz"
            )

            samtools_sort_cl = self.get_command_line_from_mock(run_core_mock, 10)
            assert (
                samtools_sort_cl
                == f"samtools sort -@ 8 -o {td}/custom_tmp/custom_name.subsampled.sorted.sam {td}/custom_tmp/custom_name.subsampled.sam"
            )

            bcftools_mpileup_cl = self.get_command_line_from_mock(run_core_mock, 11)
            assert (
                bcftools_mpileup_cl
                == f"bcftools mpileup -f {H37RV_genome} --threads 8 -x -I -Q 13 "
                f"-a INFO/SCR,FORMAT/SP,INFO/ADR,INFO/ADF -h100 -M10000 -o {td}/custom_tmp/custom_name.pileup.vcf "
                f"{td}/custom_tmp/custom_name.subsampled.sorted.sam"
            )

            bcftools_call_cl = self.get_command_line_from_mock(run_core_mock, 12)
            assert (
                bcftools_call_cl
                == f"bcftools call --threads 8 --ploidy 1 -V indels -m -o {td}/custom_tmp/custom_name.snps.vcf "
                f"{td}/custom_tmp/custom_name.pileup.vcf"
            )

            filter_vcf_cl = self.get_command_line_from_mock(run_core_mock, 13)
            assert (
                filter_vcf_cl
                == f"{sys.executable} {EXTERNAL_SCRIPTS_DIR}/apply_filters.py -P --verbose --overwrite -d 5 -D 0 "
                f"-q 25 -s 1 -b 0 -m 0 -r 0 -V 1e-05 -G 0 -K 0.9 -M 30 -x 0.2 "
                f"-o {td}/custom_name.snps.filtered.bcf -i {td}/custom_tmp/custom_name.snps.vcf"
            )

            generate_consensus_cl = self.get_command_line_from_mock(run_core_mock, 14)
            assert (
                generate_consensus_cl
                == f"{sys.executable} {EXTERNAL_SCRIPTS_DIR}/consensus.py --sample-id custom_name --verbose --ignore all "
                f"--het-default none -o {td}/custom_name.consensus.fa -i {td}/custom_name.snps.filtered.bcf "
                f"-f {H37RV_genome} -m {H37RV_mask}"
            )

    @patch.object(
        ExternalTool,
        ExternalTool._run_core.__name__,
        side_effect=["", "", subprocess.CalledProcessError(1, "minimap2")],
    )
    def test_partial_execution___minimum_params___index_sorted_decontaminated_bam_fails___checks_fail_happens_and_previous_tools_called_correctly(
        self, run_core_mock, ensure_decontamination_db_is_available_mock, tmp_path
    ):
        runner = CliRunner()
        with runner.isolated_filesystem(temp_dir=tmp_path) as td:
            td = Path(td)
            infile = td / "in.fq"
            with open(infile, "w") as fp:
                fp.write("@r1\nACGT\n+$$$%\n")
            opts = [
                "process",
                str(infile),
                "--cleanup",
                "-o",
                str(td),
            ]  # asks for cleanup, but it should not cleanup as this is a run that fails
            result = runner.invoke(main_cli, opts)

            # check if tbpore indeed failed
            assert result.exit_code == 1
            index_sorted_decontaminated_bam_cl = self.get_command_line_from_mock(
                run_core_mock, 2
            )
            assert (
                b"Error calling "
                + index_sorted_decontaminated_bam_cl.encode("utf-8")
                + b" (return code 1)"
                in result.stdout_bytes
            )

            # check if all tools until minimap2 were called correctly
            assert run_core_mock.call_count == 3

            map_decontamination_db_cl = self.get_command_line_from_mock(
                run_core_mock, 0
            )
            assert (
                map_decontamination_db_cl
                == f"minimap2 -aL2 -x map-ont -t 1 -o {td}/{TMP_NAME}/in.decontaminated.sam {DECONTAMINATION_DB_INDEX} {td}/{TMP_NAME}/in.fq.gz"
            )

            sort_decontaminated_sam_cl = self.get_command_line_from_mock(
                run_core_mock, 1
            )
            assert (
                sort_decontaminated_sam_cl
                == f"samtools sort -@ 1 -o {td}/{TMP_NAME}/in.decontaminated.sorted.bam {td}/{TMP_NAME}/in.decontaminated.sam"
            )

            assert (
                index_sorted_decontaminated_bam_cl
                == f"samtools index -@ 1 {td}/{TMP_NAME}/in.decontaminated.sorted.bam"
            )

            # check if tmp not removed
            tbpore_tmp = td / TMP_NAME
            assert tbpore_tmp.exists()


@patch.object(ExternalTool, ExternalTool._run_core.__name__)
@patch("tbpore.tbpore.ensure_decontamination_db_is_available")
class TestCLICleanup:
    def test_no_cleanup(
        self, ensure_decontamination_db_is_available_mock, run_core_mock, tmp_path
    ):
        sample = "sam"
        opts = ["process", "--no-cleanup", "-S", sample]
        runner = CliRunner()
        with runner.isolated_filesystem(temp_dir=tmp_path) as td:
            td = Path(td)
            infile = td / "in.fq"
            with open(infile, "w") as fp:
                fp.write("@r1\nACGT\n+$$$%\n")

            opts.extend(["-o", str(td)])
            opts.extend([str(infile)])
            result = runner.invoke(main_cli, opts)
            assert result.exit_code == 0
            tbpore_tmp = td / TMP_NAME
            assert tbpore_tmp.exists()

    def test_with_cleanup(
        self, ensure_decontamination_db_is_available_mock, run_core_mock, tmp_path
    ):
        sample = "sam"
        opts = ["process", "--cleanup", "-S", sample]
        runner = CliRunner()
        with runner.isolated_filesystem(temp_dir=tmp_path) as td:
            td = Path(td)
            infile = td / "in.fq"
            with open(infile, "w") as fp:
                fp.write("@r1\nACGT\n+$$$%\n")

            opts.extend(["-o", str(td)])
            opts.extend([str(infile)])
            result = runner.invoke(main_cli, opts)
            assert result.exit_code == 0
            tbpore_tmp = td / TMP_NAME
            assert not tbpore_tmp.exists()


@patch.object(ExternalTool, ExternalTool._run_core.__name__)
@patch("tbpore.tbpore.ensure_decontamination_db_is_available")
class TestInputConcatenation:
    def test_no_input___fails(
        self, ensure_decontamination_db_is_available_mock, run_core_mock, tmp_path
    ):
        opts = ["process"]
        runner = CliRunner()
        with runner.isolated_filesystem(temp_dir=tmp_path):
            result = runner.invoke(main_cli, opts)
            assert result.exit_code == 2
            assert b"No INPUT files given" in result.stdout_bytes

    def test_single_file(
        self, ensure_decontamination_db_is_available_mock, run_core_mock, tmp_path
    ):
        sample = "sam"
        opts = ["process", "-D", "-S", sample]
        runner = CliRunner()
        expected_fq = "@r1\nACGT\n+$$$%\n"
        with runner.isolated_filesystem(temp_dir=tmp_path) as td:
            td = Path(td)
            infile = td / "in.fq"
            with open(infile, "w") as fp:
                fp.write(expected_fq)

            opts.extend(["-o", str(td)])
            opts.extend([str(infile)])
            result = runner.invoke(main_cli, opts)
            assert result.exit_code == 0
            tbpore_concat = td / TMP_NAME / f"{sample}.fq.gz"
            assert tbpore_concat.exists()
            with gzip.open(tbpore_concat, mode="rt") as fp:
                actual = fp.read()
            assert actual == expected_fq

    def test_single_gz_file(
        self, ensure_decontamination_db_is_available_mock, run_core_mock, tmp_path
    ):
        sample = "sam"
        opts = ["process", "-D", "-S", sample]
        runner = CliRunner()
        expected_fq = "@r1\nACGT\n+$$$%\n"
        with runner.isolated_filesystem(temp_dir=tmp_path) as td:
            td = Path(td)
            infile = td / "in.fastq.gz"
            with gzip.open(infile, "wb") as fp:
                fp.write(expected_fq.encode())

            opts.extend(["-o", str(td)])
            opts.extend([str(infile)])
            result = runner.invoke(main_cli, opts)
            assert result.exit_code == 0
            tbpore_concat = td / TMP_NAME / f"{sample}.fq.gz"
            assert tbpore_concat.exists()
            with gzip.open(tbpore_concat, mode="rt") as fp:
                actual = fp.read()
            assert actual == expected_fq

    def test_multiple_files(
        self, ensure_decontamination_db_is_available_mock, run_core_mock, tmp_path
    ):
        sample = "sam"
        opts = ["process", "-D", "-S", sample]
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
            result = runner.invoke(main_cli, opts)
            assert result.exit_code == 0
            tbpore_concat = td / TMP_NAME / f"{sample}.fq.gz"
            assert tbpore_concat.exists()
            with gzip.open(tbpore_concat, mode="rt") as fp:
                actual = fp.read()

            expected = expected_fq1 + expected_fq2
            assert sorted(actual) == sorted(expected)

    def test_input_is_empty_dir___error_out(
        self, ensure_decontamination_db_is_available_mock, run_core_mock, tmp_path
    ):
        sample = "sam"
        opts = ["process", "-D", "-S", sample]
        runner = CliRunner()
        with runner.isolated_filesystem(temp_dir=tmp_path) as td:
            td = Path(td)
            opts.extend(["-o", str(td)])
            opts.extend([str(td)])
            result = runner.invoke(main_cli, opts)
            assert result.exit_code == 2
            assert (
                b"No fastq files found for the given inputs, please check your input files/dirs."
                in result.stdout_bytes
            )

    def test_input_is_dir(
        self, ensure_decontamination_db_is_available_mock, run_core_mock, tmp_path
    ):
        sample = "sam"
        opts = ["process", "-D", "-S", sample]
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
            result = runner.invoke(main_cli, opts)
            assert result.exit_code == 0
            tbpore_concat = td / TMP_NAME / f"{sample}.fq.gz"
            assert tbpore_concat.exists()
            with gzip.open(tbpore_concat, mode="rt") as fp:
                actual = fp.read()

            expected = expected_fq1 + expected_fq2
            assert sorted(actual) == sorted(expected)

    def test_input_is_dir_and_one_file_has_bad_suffix(
        self, ensure_decontamination_db_is_available_mock, run_core_mock, tmp_path
    ):
        sample = "sam"
        opts = ["process", "-D", "-S", sample]
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
            result = runner.invoke(main_cli, opts)
            assert result.exit_code == 0
            tbpore_concat = td / TMP_NAME / f"{sample}.fq.gz"
            assert tbpore_concat.exists()
            with gzip.open(tbpore_concat, mode="rt") as fp:
                actual = fp.read()

            expected = expected_fq2
            assert sorted(actual) == sorted(expected)

    def test_input_is_dir_and_files_in_dir___ensure_duplication_does_not_happen(
        self, ensure_decontamination_db_is_available_mock, run_core_mock, tmp_path
    ):
        sample = "sam"
        opts = ["process", "-D", "-S", sample]
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
            result = runner.invoke(main_cli, opts)
            assert result.exit_code == 0
            assert b"Found 2 fastq files. Joining them..." in result.stdout_bytes

            tbpore_concat = td / TMP_NAME / f"{sample}.fq.gz"
            assert tbpore_concat.exists()
            with gzip.open(tbpore_concat, mode="rt") as fp:
                actual = fp.read()

            expected = expected_fq1 + expected_fq2
            assert sorted(actual) == sorted(expected)
