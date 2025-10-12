#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import textwrap
import statistics
import random
from random import randint
from pathlib import Path
from typing import Iterator, Dict, List, NamedTuple
import matplotlib
import matplotlib.pyplot as plt
import networkx as nx
from networkx import (
    DiGraph,
    all_simple_paths,
    lowest_common_ancestor,
    has_path,
)

random.seed(9001)

matplotlib.use("Agg")

__author__ = "Anaïs DELASSUS"
__copyright__ = "Universite Paris Cité"
__credits__ = ["Anaïs DELASSUS"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Anaïs DELASSUS"
__email__ = "anais.delassus@etu.u-paris.fr"
__status__ = "Developpement"


def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def get_arguments():  # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(
        description=__doc__, usage=f"{sys.argv[0]} -h"
    )
    parser.add_argument(
        "-i", dest="fastq_file", type=isfile, required=True, help="Fastq file"
    )
    parser.add_argument(
        "-k", dest="kmer_size", type=int, default=22,
        help="k-mer size (default 22)"
    )
    parser.add_argument(
        "-o",
        dest="output_file",
        type=Path,
        default=Path(os.curdir + os.sep + "contigs.fasta"),
        help="Output contigs in fasta file (default contigs.fasta)",
    )
    parser.add_argument(
        "-f", dest="graphimg_file", type=Path,
        help="Save graph as an image (png)"
    )
    return parser.parse_args()


def read_fastq(fastq_file: Path) -> Iterator[str]:
    """Extract reads from fastq files.

    :param fastq_file: (Path) Path to the fastq file.
    :return: A generator object that iterate the read sequences.
    """
    with open(fastq_file, "r", encoding="utf-8") as fq:
        for _, seq, _, _ in zip(*[fq]*4):
            yield seq.strip()


def cut_kmer(read: str, kmer_size: int) -> Iterator[str]:
    """Cut read into kmers of size kmer_size.

    :param read: (str) Sequence of a read.
    :return: A generator object that provides the kmers (str) of size 
    kmer_size.
    """
    for i in range(len(read) - kmer_size + 1):
        yield read[i: i + kmer_size]


def build_kmer_dict(fastq_file: Path, kmer_size: int) -> Dict[str, int]:
    """Build a dictionnary object of all kmer occurrences in the fastq file

    :param fastq_file: (str) Path to the fastq file.
    :return: A dictionnary object that identify all kmer occurrences.
    """
    occ_kmer = {}
    for read in read_fastq(fastq_file):
        for kmer in cut_kmer(read, kmer_size):
            if kmer in occ_kmer:
                occ_kmer[kmer] += 1
            else:
                occ_kmer[kmer] = 1
    return occ_kmer


def build_graph(kmer_dict: Dict[str, int]) -> DiGraph:
    """Build the debruijn graph

    :param kmer_dict: A dictionnary object that identify all kmer occurrences.
    :return: A directed graph (nx) of all kmer substring and weight (occ).
    """
    graph = DiGraph()
    for kmer, occ in kmer_dict.items():
        prefixe = kmer[:-1]
        suffixe = kmer[1:]
        graph.add_edge(prefixe, suffixe, weight=occ)
    return graph


def remove_paths(
    graph: DiGraph,
    path_list: List[List[str]],
    delete_entry_node: bool,
    delete_sink_node: bool,
) -> DiGraph:
    """Remove a list of path in a graph. A path is set of connected node in
    the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param delete_entry_node: (boolean) True->We remove the 1st node of a path
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    nodes_to_remove = set()
    for path in path_list:
        if not path:
            continue
        if len(path) == 1:
            if delete_entry_node and delete_sink_node:
                nodes_to_remove.add(path[0])
            continue
        start_idx = 0
        end_idx = len(path) - 1
        del_start = start_idx if delete_entry_node else start_idx + 1
        del_end = end_idx if delete_sink_node else end_idx - 1
        for i in range(del_start, del_end + 1):
            if 0 <= i < len(path):
                nodes_to_remove.add(path[i])
    graph.remove_nodes_from(nodes_to_remove)
    return graph


class PathOptions(NamedTuple):
    """Evite de mettre trop d'argument dans select_best_path()."""
    delete_entry_node: bool = False
    delete_sink_node: bool = False


def choose_best_index(lengths: List[int], weights: List[float]) -> int:
    """Return the index of the best path based on weight or length"""
    n = len(weights)
    try:
        w_stdev = statistics.stdev(weights) if n >= 2 else 0.0
    except statistics.StatisticsError:
        w_stdev = 0.0
    if w_stdev > 0:
        return max(range(n), key=lambda i: weights[i])
    try:
        l_stdev = statistics.stdev(lengths) if n >= 2 else 0.0
    except statistics.StatisticsError:
        l_stdev = 0.0
    if l_stdev > 0:
        return max(range(n), key=lambda i: lengths[i])
    return randint(0, n - 1)


def select_best_path(
    graph: DiGraph,
    path_list: List[List[str]],
    path_length: List[int],
    weight_avg_list: List[float],
    options: PathOptions = PathOptions()
) -> DiGraph:
    """Select the best path between different paths"""
    delete_entry_node = options.delete_entry_node
    delete_sink_node = options.delete_sink_node
    if not path_list or len(path_list) <= 1:
        return graph
    if len(path_length) != len(path_list) or\
            len(weight_avg_list) != len(path_list):
        raise ValueError("path_length and weight_avg_list must match path_list\
            length")
    best_index = choose_best_index(path_length, weight_avg_list)
    bad_paths = [p for idx, p in enumerate(path_list) if idx != best_index]
    if bad_paths:
        graph = remove_paths(graph, bad_paths, delete_entry_node,
                             delete_sink_node)
    return graph


def path_average_weight(graph: DiGraph, path: List[str]) -> float:
    """Compute the weight of a path

    :param graph: (nx.DiGraph) A directed graph object
    :param path: (list) A path consist of a list of nodes
    :return: (float) The average weight of a path
    """
    return statistics.mean(
        [d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)]
    )


def solve_bubble(graph: DiGraph, ancestor_node: str, descendant_node: str
                 ) -> DiGraph:
    """Explore and solve bubble issue

    :param graph: (nx.DiGraph) A directed graph object
    :param ancestor_node: (str) An upstream node in the graph
    :param descendant_node: (str) A downstream node in the graph
    :return: (nx.DiGraph) A directed graph object
    """
    path_list = list(nx.all_simple_paths(
        graph, source=ancestor_node, target=descendant_node))
    if len(path_list) <= 1:
        return graph
    path_length_list: List[int] = []
    weight_avg_list: List[float] = []
    for path in path_list:
        path_length_list.append(len(path))
        edges_data = list(graph.subgraph(path).edges(data=True))
        weights = [d["weight"] for (_, _, d) in edges_data if "weight" in d]
        weight_avg = statistics.mean(weights) if weights else 0
        weight_avg_list.append(weight_avg)
    graph = select_best_path(
        graph,
        path_list,
        path_length_list,
        weight_avg_list
    )
    return graph


def simplify_bubbles(graph: DiGraph) -> DiGraph:
    """
    Detect and explode bubbles in a directed graph.
    :param graph: (nx.DiGraph) A directed graph
    :return: (nx.DiGraph) Graph with bubbles simplified
    """
    bubble_found = False
    bubble_start = None
    bubble_end = None
    for node in graph.nodes:
        predecessors = list(graph.predecessors(node))
        if len(predecessors) > 1:
            for i, pred_i in enumerate(predecessors):
                for pred_j in predecessors[i + 1:]:
                    anc = lowest_common_ancestor(graph, pred_i, pred_j)
                    if anc is not None:
                        bubble_found = True
                        bubble_start = anc
                        bubble_end = node
                        break
                if bubble_found:
                    break
        if bubble_found:
            break
    if bubble_found:
        graph = solve_bubble(graph, bubble_start, bubble_end)
        graph = simplify_bubbles(graph)
    return graph


def solve_entry_tips(graph: DiGraph, starting_nodes: List[str]) -> DiGraph:
    """Remove entry tips (short dead-end branches at graph entry)"""
    tip_found = False
    for node in starting_nodes:
        successors = list(graph.successors(node))
        if len(successors) > 1:
            # Plusieurs branches partent du même noeud d'entrée
            path_list = []
            path_length_list = []
            weight_avg_list = []
            for succ in successors:
                paths = list(all_simple_paths(graph, source=node, target=succ))
                for path in paths:
                    path_list.append(path)
                    path_length_list.append(len(path))
                    weight_avg_list.append(path_average_weight(graph, path))
            if path_list:
                options = PathOptions(delete_entry_node=True,
                                      delete_sink_node=False)
                graph = select_best_path(
                    graph, path_list, path_length_list, weight_avg_list,
                    options=options)
                tip_found = True
    return graph if tip_found else graph


def solve_out_tips(graph: DiGraph, ending_nodes: List[str]) -> DiGraph:
    """Remove out tips (short dead-end branches at graph exit)"""
    tip_found = False
    for node in ending_nodes:
        predecessors = list(graph.predecessors(node))
        if len(predecessors) > 1:
            # Plusieurs branches arrivent au même noeud de sortie
            path_list = []
            path_length_list = []
            weight_avg_list = []
            for pred in predecessors:
                paths = list(all_simple_paths(graph, source=pred, target=node))
                for path in paths:
                    path_list.append(path)
                    path_length_list.append(len(path))
                    weight_avg_list.append(path_average_weight(graph, path))
            if path_list:
                options = PathOptions(delete_entry_node=False,
                                      delete_sink_node=True)
                graph = select_best_path(
                    graph, path_list, path_length_list, weight_avg_list,
                    options=options
                )
                tip_found = True
    return graph if tip_found else graph


def get_starting_nodes(graph: DiGraph) -> List[str]:
    """Get nodes without predecessors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without predecessors
    """
    noeuds_entree = []
    for node in graph.nodes():
        if len(list(graph.predecessors(node))) == 0:
            noeuds_entree.append(node)
    return noeuds_entree


def get_sink_nodes(graph: DiGraph) -> List[str]:
    """Get nodes without successors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without successors
    """
    noeuds_sortie = []
    for node in graph.nodes():
        if len(list(graph.successors(node))) == 0:
            noeuds_sortie.append(node)
    return noeuds_sortie


def get_contigs(
    graph: DiGraph, starting_nodes: List[str], ending_nodes: List[str]
) -> List:
    """Extract the contigs from the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param starting_nodes: (list) A list of nodes without predecessors
    :param ending_nodes: (list) A list of nodes without successors
    :return: (list) List of [contiguous sequence and their length]
    """
    contigs = []
    for entree in starting_nodes:
        for sortie in ending_nodes:
            if has_path(graph, entree, sortie):
                for path in all_simple_paths(graph, entree, sortie):
                    contig = path[0]
                    for node in path[1:]:
                        contig += node[-1]
                    contigs.append((contig, len(contig)))
    return contigs


def save_contigs(contigs_list: List[str], output_file: Path) -> None:
    """Write all contigs in fasta format

    :param contig_list: (list) List of [contiguous sequence and their length]
    :param output_file: (Path) Path to the output file
    """
    with open(output_file, "w", encoding="utf-8") as output:
        for i, (contig, taille) in enumerate(contigs_list):
            output.write(f">contig_{i} len={taille}\n")
            output.write(textwrap.fill(contig, width=80) + "\n")


def draw_graph(graph: DiGraph, graphimg_file: Path) -> None:
    """Draw the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param graphimg_file: (Path) Path to the output file
    """
    plt.figure(figsize=(12, 12))
    # spring_layout donne une disposition plus lisible
    pos = nx.spring_layout(graph, k=0.5, iterations=100)

    # séparer arêtes fortes/faibles selon poids
    elarge = [(u, v) for u, v, d in graph.edges(
        data=True) if d.get("weight", 0) > 3]
    esmall = [(u, v) for u, v, d in graph.edges(
        data=True) if d.get("weight", 0) <= 3]

    nx.draw_networkx_nodes(graph, pos, node_size=10, alpha=0.7)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall,
                           width=0.5, alpha=0.2, edge_color="gray")
    nx.draw_networkx_edges(graph, pos, edgelist=elarge,
                           width=1.5, alpha=0.8, edge_color="red")

    plt.axis("off")
    plt.tight_layout()
    plt.savefig(graphimg_file, dpi=300)
    plt.close()


# ==============================================================
# Main program
def main() -> None:
    """Main program function (lecture, nettoyage du graphe, 
    extraction des contigs)."""

    # 1. Lecture du fichier et construction du graphe
    args = get_arguments()
    kmer_dict = build_kmer_dict(args.fastq_file, args.kmer_size)
    print(f"{len(kmer_dict)} k-mers générés.")
    graph = build_graph(kmer_dict)
    print(f"Graphe construit avec {graph.number_of_nodes()} noeuds et {
          graph.number_of_edges()} arêtes.")

    # 2. Résolution des bulles
    graph = simplify_bubbles(graph)
    print(f"Graphe après simplification des bulles : {
          graph.number_of_nodes()} noeuds, {graph.number_of_edges()} arêtes.")

    # Identifier les noeuds d'entrée et de sortie
    starting_nodes = get_starting_nodes(graph)
    ending_nodes = get_sink_nodes(graph)
    print(f"{len(starting_nodes)} noeuds d'entrée, {
          len(ending_nodes)} noeuds de sortie.")

    # 3. Résolution des pointes d’entrée et de sortie
    if starting_nodes:
        graph = solve_entry_tips(graph, starting_nodes)
    if ending_nodes:
        graph = solve_out_tips(graph, ending_nodes)

    # Recalculer les noeuds après modification du graphe
    starting_nodes = get_starting_nodes(graph)
    ending_nodes = get_sink_nodes(graph)
    print(f"Après nettoyage des tips : {graph.number_of_nodes()} noeuds, {
          graph.number_of_edges()} arêtes.")
    print(f"{len(starting_nodes)} noeuds d'entrée, {
          len(ending_nodes)} noeuds de sortie.")

    # 4. Ecriture du/des contigs
    contigs = get_contigs(graph, starting_nodes, ending_nodes)
    print(f"{len(contigs)} contigs générés.")
    save_contigs(contigs, args.output_file)
    print(f"Contigs sauvegardés dans {args.output_file}.")

    # Optionnel : sauvegarde du graphe en image
    if args.graphimg_file:
        draw_graph(graph, args.graphimg_file)
        print(f"Graphe sauvegardé dans {args.graphimg_file}.")


if __name__ == "__main__":  # pragma: no cover
    main()
