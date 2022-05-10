from pathlib import Path
import subprocess

# configs
samples = ["mada_1-19", "mada_1-8", "mada_1-6", "mada_103", "mada_1-1", "mada_130", "mada_132", "mada_128", "mada_2-1", "mada_1-51", "mada_1-20", "mada_2-31", "mada_109", "mada_112", "mada_1-5", "mada_107", "mada_1-40", "mada_111", "mada_1-14", "mada_1-46", "mada_136", "mada_1-54", "mada_1-25", "mada_118", "mada_129", "mada_1-18", "mada_151", "mada_134", "mada_1-3", "mada_1-44", "mada_1-15", "mada_1-47", "mada_2-53", "mada_1-16", "mada_1-2", "mada_154", "mada_104", "mada_115", "mada_126", "mada_1-33", "mada_102", "mada_127", "mada_125", "mada_1-43", "mada_137", "mada_117", "mada_131", "mada_2-46", "mada_1-17", "mada_124", "mada_121", "mada_141", "mada_2-42", "mada_1-28", "mada_152", "mada_1-48", "mada_2-34", "mada_123", "mada_106", "mada_140", "mada_1-30", "mada_1-50", "mada_139", "mada_1-41", "mada_2-25", "mada_105", "mada_144", "mada_1-32", "mada_1-22", "mada_1-11", "mada_1-10", "mada_110", "mada_120", "mada_113", "mada_1-38", "mada_122", "mada_1-7", "mada_116", "mada_142", "mada_133", "mada_148", "mada_1-13", "mada_1-53", "mada_135", "mada_1-12", "mada_2-50", "mada_1-21", "mada_1-39", "mada_1-36", "mada_143", "mada_150"]
H2H_path = Path("/hps/nobackup/iqbal/mbhall/tech_wars/data/QC/filtered/madagascar/nanopore")
tbpore_path = Path("/hps/nobackup/iqbal/leandro/tbpore2/tbpore/pipelines/snakemake/output_human_decon")

sample_to_nb_of_diffs = {}
for sample in samples:
    h2h_keep_reads = H2H_path/sample/"keep.reads"
    tbpore_keep_reads = tbpore_path/sample/".tbpore"/f"{sample}.decontaminated.filter/keep.reads"
    diff_out = subprocess.check_output(f"diff --suppress-common-lines -y <(sort {h2h_keep_reads}) <(sort {tbpore_keep_reads})  | wc -l", shell=True)
    nb_of_diffs = int(diff_out.strip())
    sample_to_nb_of_diffs[sample] = nb_of_diffs

print("Differences:")
for sample, nb_of_diffs in sample_to_nb_of_diffs.items():
    print(f"{sample} {nb_of_diffs}")
print(f"Total {sum(sample_to_nb_of_diffs.values())}")
