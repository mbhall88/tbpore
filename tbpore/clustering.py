"""
Module responsible for clustering tbpore process consensus sequences.
Totally based on https://github.com/mbhall88/head_to_head_pipeline/blob/bcbc84971342a26cd0a9f0ad8df4f01dcf35c01c/analysis/transmission_clustering/eda/clustering.ipynb
"""

from itertools import chain
from pathlib import Path
from typing import Generator, List, Set

import networkx as nx
import numpy as np
import pandas as pd

DELIM = ","
PAIR_IDX = ("sample1", "sample2")


class AsymmetrixMatrixError(Exception):
    pass


def load_matrix(fpath, delim: str = DELIM, name: str = "") -> pd.Series:
    matrix = []
    with open(fpath) as instream:
        header = next(instream).rstrip()
        names = np.array(header.split(delim)[1:])
        idx = np.argsort(names)
        sorted_names = names[idx]
        for row in map(str.rstrip, instream):
            # sort row according to the name sorting
            sorted_row = np.array(row.split(delim)[1:], dtype=int)[idx]
            matrix.append(sorted_row)
    sorted_matrix = np.array(matrix)[idx]
    n_samples = len(sorted_names)
    diagonal_is_zero = all(sorted_matrix[i, i] == 0 for i in range(n_samples))
    if not diagonal_is_zero:
        raise AsymmetrixMatrixError("Distance matrix diagonal is not all zero")

    matrix_is_symmetric = np.allclose(sorted_matrix, sorted_matrix.T)
    if not matrix_is_symmetric:
        raise AsymmetrixMatrixError("Distance matrix is not symmetric")

    mx = pd.DataFrame(sorted_matrix, columns=sorted_names, index=sorted_names)
    # remove the lower triangle of the matrix and the middle diagonal
    mx = mx.where(np.triu(np.ones(mx.shape), k=1).astype(bool))
    mx = mx.stack().rename(name).astype(int)
    mx = mx.rename_axis(PAIR_IDX)

    return mx


def matrix_to_graph(
    mx: pd.Series, threshold: int, include_singletons: bool = False
) -> nx.Graph:
    edges = [(s1, s2, dist) for (s1, s2), dist in mx.items() if dist <= threshold]
    graph = nx.Graph()
    graph.add_weighted_edges_from(edges)
    if include_singletons:
        samples = set()
        for u in chain.from_iterable(mx.index):
            if u not in samples:
                graph.add_node(u)
                samples.add(u)
    return graph


def sort_clusters(clusters: Generator[Set[str], None, None]) -> List[List[str]]:
    # gets a list of sorted clusters
    clusters = [sorted(list(cluster)) for cluster in clusters]

    # sort by size and then IDs (in case of draws)
    clusters = sorted(clusters, key=lambda cluster: (-len(cluster), str(cluster)))

    return clusters


def get_clusters(psdm_matrix: Path, clustering_threshold: int) -> List[List[str]]:
    ont_mtx = load_matrix(psdm_matrix, name="nanopore")
    ont_graph = matrix_to_graph(
        ont_mtx, threshold=clustering_threshold, include_singletons=True
    )

    clusters = nx.connected_components(ont_graph)

    # sort clusters to guarantee determinism
    clusters = sort_clusters(clusters)
    return clusters


def get_formatted_clusters(clusters: List[List[str]]) -> str:
    clusters_as_strs = []
    for cluster_index, cluster in enumerate(clusters):
        cluster_as_str = "\t".join(cluster)
        cluster_description = f"Cluster #{cluster_index+1}:\t{cluster_as_str}"
        clusters_as_strs.append(cluster_description)
    return "\n".join(clusters_as_strs)


def produce_clusters(psdm_matrix: Path, threshold: int, outdir: Path) -> None:
    clusters = get_clusters(psdm_matrix, threshold)
    clusters_file = outdir / "clusters.txt"
    with open(clusters_file, "w") as clusters_fh:
        print(get_formatted_clusters(clusters), file=clusters_fh)
