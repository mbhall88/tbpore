from glob import glob
from pathlib import Path
import subprocess


def get_samples(glob_pattern):
    consensus_files = (Path(file).resolve() for file in glob(glob_pattern))
    return {file.name.replace(".consensus.fa", ""): file for file in consensus_files}


def compare(h2h_glob_pattern, tbpore_glob_pattern, output):
    h2h_samples = get_samples(h2h_glob_pattern)
    tbpore_samples = get_samples(tbpore_glob_pattern)
    common_samples = set(h2h_samples.keys()).intersection(set(tbpore_samples.keys()))

    for sample in common_samples:
        h2h_consensus = h2h_samples[sample]
        tbpore_consensus = tbpore_samples[sample]
        subprocess.check_call(f"psdm -l -s -i -P -t 1 {h2h_consensus} {tbpore_consensus} >>{output}  2>/dev/null", shell=True)


def __main__():
    with open(snakemake.output.consensus_comparison, "w"):
        pass
    compare(h2h_glob_pattern=snakemake.params.h2h_glob_pattern,
            tbpore_glob_pattern=f"{snakemake.params.tbpore_output}/*/*.consensus.fa",
            output=snakemake.output.consensus_comparison)


__main__()
