import json
from glob import glob
from pathlib import Path

import dictdiffer


def get_samples(glob_pattern):
    mykrobe_files = (Path(file).resolve() for file in glob(glob_pattern))
    return {file.name.replace(".mykrobe.json", ""): file for file in mykrobe_files}


def cleanup_susceptiblity_dict(dict):
    for drug in dict:
        if "called_by" in dict[drug]:
            del dict[drug]["called_by"]
    for drug in dict:
        assert "called_by" not in dict[drug]


def compare(h2h_mykrobe_glob_pattern, tbpore_glob_pattern, output_fh, output_simple_fh):
    h2h_samples = get_samples(h2h_mykrobe_glob_pattern)
    tbpore_samples = get_samples(tbpore_glob_pattern)
    common_samples = set(h2h_samples.keys()).intersection(set(tbpore_samples.keys()))
    for sample in common_samples:
        h2h_mykrobe = h2h_samples[sample]
        tbpore_mykrobe = tbpore_samples[sample]

        with open(h2h_mykrobe) as h2h_mykrobe_fh, open(
            tbpore_mykrobe
        ) as tbpore_mykrobe_fh:
            h2h_mykrobe_json = json.load(h2h_mykrobe_fh)
            tbpore_mykrobe_json = json.load(tbpore_mykrobe_fh)

        print(f"Diffing mykrobe jsons for {sample}", file=output_fh)
        h2h_susceptibility = h2h_mykrobe_json[sample]["susceptibility"]
        tbpore_susceptibility = tbpore_mykrobe_json[sample]["susceptibility"]
        for diff in list(dictdiffer.diff(h2h_susceptibility, tbpore_susceptibility)):
            print(diff, file=output_fh)

        print(f"Diffing mykrobe jsons for {sample}", file=output_simple_fh)
        cleanup_susceptiblity_dict(h2h_susceptibility)
        cleanup_susceptiblity_dict(tbpore_susceptibility)
        for diff in list(dictdiffer.diff(h2h_susceptibility, tbpore_susceptibility)):
            print(diff, file=output_simple_fh)


with open(snakemake.output.mykrobe_comparison, "w") as output_fh, open(
    snakemake.output.mykrobe_comparison_simple, "w"
) as output_simple_fh:
    compare(
        h2h_mykrobe_glob_pattern=snakemake.params.h2h_mykrobe_glob_pattern,
        tbpore_glob_pattern=f"{snakemake.params.tbpore_output}/*/*.mykrobe.json",
        output_fh=output_fh,
        output_simple_fh=output_simple_fh,
    )
