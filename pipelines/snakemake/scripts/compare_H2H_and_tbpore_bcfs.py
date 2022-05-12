from glob import glob
from pathlib import Path

from cyvcf2 import VCF


def get_samples(glob_pattern):
    bcf_files = (Path(file).resolve() for file in glob(glob_pattern))
    return {file.name.replace(".snps.filtered.bcf", ""): file for file in bcf_files}


def get_var_with_no_qual(variant):
    variant = str(variant)
    split_variant = variant.split("\t")
    split_variant_with_no_qual = split_variant[:5] + split_variant[6:]
    return "\t".join(split_variant_with_no_qual)


def compare(h2h_bcf_glob_pattern, tbpore_glob_pattern, output):
    h2h_samples = get_samples(h2h_bcf_glob_pattern)
    tbpore_samples = get_samples(tbpore_glob_pattern)
    common_samples = set(h2h_samples.keys()).intersection(set(tbpore_samples.keys()))

    with open(output, "w") as output_fh:
        for sample in common_samples:
            h2h_bcf = VCF(h2h_samples[sample])
            tbpore_bcf = VCF(tbpore_samples[sample])
            nb_of_different_variants = 0
            for h2h_variant, tbpore_variant in zip(h2h_bcf, tbpore_bcf):
                h2h_variant_with_no_qual = get_var_with_no_qual(h2h_variant)
                tbpore_variant_with_no_qual = get_var_with_no_qual(tbpore_variant)
                variants_are_different = (
                    h2h_variant_with_no_qual != tbpore_variant_with_no_qual
                )
                if variants_are_different:
                    print(
                        f"Different VCF records in sample {sample}:\n{h2h_variant}{tbpore_variant}",
                        file=output_fh,
                    )
                    nb_of_different_variants += 1
            print(
                f"Sample {sample}: {nb_of_different_variants} different variants",
                file=output_fh,
            )


def __main__():
    compare(
        h2h_bcf_glob_pattern=snakemake.params.h2h_bcf_glob_pattern,
        tbpore_glob_pattern=f"{snakemake.params.tbpore_output}/*/*.snps.filtered.bcf",
        output=snakemake.output.bcfs_comparison,
    )


__main__()
