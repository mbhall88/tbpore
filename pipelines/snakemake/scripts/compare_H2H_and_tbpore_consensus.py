from glob import glob
from pathlib import Path
import gzip
from Bio import SeqIO


def get_samples(glob_pattern):
    consensus_files = (Path(file).resolve() for file in glob(glob_pattern))
    return {file.name.replace(".consensus.fa", ""): file for file in consensus_files}


def base_to_key(base, ref_base):
    is_ref = base == ref_base
    is_N = not is_ref and base == "N"
    is_alt = not is_N
    if is_ref:
        return "ref"
    if is_N:
        return "N"
    if is_alt:
        return "alt"


def compare_core(h2h_consensus, tbpore_consensus, H37RV_genome):
    with open(h2h_consensus) as h2h_consensus_fh, \
         open(tbpore_consensus) as tbpore_consensus_fh, \
         gzip.open(H37RV_genome) as H37RV_genome_fh:
        h2h_consensus_records = list(SeqIO.parse(h2h_consensus_fh, "fasta"))
        tbpore_consensus_records = list(SeqIO.parse(tbpore_consensus_fh, "fasta"))
        H37RV_genome_records = list(SeqIO.parse(H37RV_genome_fh, "fasta"))
        assert {len(h2h_consensus_records), len(tbpore_consensus_records), len(H37RV_genome_records)} == {1}

    h2h_consensus_seq = str(h2h_consensus_records[0].seq)
    tbpore_consensus_seq = str(tbpore_consensus_records[0].seq)
    H37RV_genome_seq = str(H37RV_genome_records[0].seq)
    assert len({len(h2h_consensus_seq), len(tbpore_consensus_seq), len(H37RV_genome_seq)}) == 1

    stats_to_count = {
        "equal": 0,
        "(ref, alt)": 0,
        "(ref, N)": 0,
        "(alt, ref)": 0,
        "(alt, alt)": 0,
        "(alt, N)": 0,
        "(N, ref)": 0,
        "(N, alt)": 0,
    }
    for h2h_base, tbpore_base, ref_base in zip(h2h_consensus_seq, tbpore_consensus_seq, H37RV_genome_seq):
        if h2h_base == tbpore_base:
            key = "equal"
        else:
            h2h_key = base_to_key(h2h_base, ref_base)
            tbpore_key = base_to_key(tbpore_base, ref_base)
            key = f"({h2h_key}, {tbpore_key})"
        stats_to_count[key] += 1

    return stats_to_count


def pretty_print_stats(stats_to_count):
    stats_to_count_sorted = dict(sorted(stats_to_count.items()))
    pp_stats = []
    for stats, count in stats_to_count_sorted.items():
        pp_stats.append(f"{stats} = {count}")
    return "\n".join(pp_stats)


def compare(h2h_consensus_glob_pattern, tbpore_glob_pattern, output):
    h2h_samples = get_samples(h2h_consensus_glob_pattern)
    tbpore_samples = get_samples(tbpore_glob_pattern)
    common_samples = set(h2h_samples.keys()).intersection(set(tbpore_samples.keys()))

    with open(output, "w") as output_fh:
        for sample in common_samples:
            h2h_consensus = h2h_samples[sample]
            tbpore_consensus = tbpore_samples[sample]
            stats_to_count = compare_core(h2h_consensus, tbpore_consensus, snakemake.params.H37RV_genome)
            print(f"Comparing {sample}:\n{pretty_print_stats(stats_to_count)}", file=output_fh)


def __main__():
    with open(snakemake.output.consensus_comparison, "w"):
        pass
    compare(h2h_consensus_glob_pattern=snakemake.params.h2h_consensus_glob_pattern,
            tbpore_glob_pattern=f"{snakemake.params.tbpore_output}/*/*.consensus.fa",
            output=snakemake.output.consensus_comparison)


__main__()
