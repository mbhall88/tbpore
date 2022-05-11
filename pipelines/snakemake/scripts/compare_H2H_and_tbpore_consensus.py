import subprocess
from glob import glob
from pathlib import Path


def get_samples(glob_pattern):
    consensus_files = (Path(file).resolve() for file in glob(glob_pattern))
    return {file.name.replace(".consensus.fa", ""): file for file in consensus_files}


def compare(h2h_consensus_glob_pattern, tbpore_glob_pattern, output):
    h2h_samples = get_samples(h2h_consensus_glob_pattern)
    tbpore_samples = get_samples(tbpore_glob_pattern)
    common_samples = set(h2h_samples.keys()).intersection(set(tbpore_samples.keys()))

    for sample in common_samples:
        h2h_consensus = h2h_samples[sample]
        tbpore_consensus = tbpore_samples[sample]
        subprocess.check_call(f"diff -qs {h2h_consensus} {tbpore_consensus} >> {output}")


def __main__():
    with open(snakemake.output.consensus_comparison, "w"):
        pass
    compare(
        h2h_consensus_glob_pattern=snakemake.params.h2h_consensus_glob_pattern,
        tbpore_glob_pattern=f"{snakemake.params.tbpore_output}/*/*.consensus.fa",
        output=snakemake.output.consensus_comparison,
    )


__main__()
