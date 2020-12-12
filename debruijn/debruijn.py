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
import statistics
import random
import networkx as nx

random.seed(9001)


__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
    :Parameters:
        path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
    Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(
        description=__doc__, usage="{0} -h".format(sys.argv[0])
    )
    parser.add_argument(
        "-i", dest="fastq_file", type=isfile, required=True, help="Fastq file"
    )
    parser.add_argument(
        "-k", dest="kmer_size", type=int, default=21, help="K-mer size (default 21)"
    )
    parser.add_argument(
        "-o",
        dest="output_file",
        type=str,
        default=os.curdir + os.sep + "contigs.fasta",
        help="Output contigs in fasta file",
    )
    return parser.parse_args()


def read_fastq(fastq_file):
    """ Reads fastq file.
    """
    with open(fastq_file, "r") as file:
        for _ in file:
            yield next(file)[:-1]  # remove '\n' character
            next(file)  # skip useless line
            next(file)  # skip useless line


def cut_kmer(read, kmer_size):
    """ Returns a generator of k-mer given a sequence.
    """
    for i in range(0, len(read) - kmer_size + 1):
        yield read[i : i + kmer_size]
        if i >= len(read) - kmer_size + 1:
            return


def build_kmer_dict(fastq_file, kmer_size):
    """ Builds a dictionary with kmer as keys and their number
    of occurences as values.
    """
    kmer_dict = {}
    for sequence in read_fastq(fastq_file):
        for kmer in cut_kmer(sequence, kmer_size):
            if kmer not in kmer_dict.keys():
                kmer_dict[kmer] = 1
            else:
                kmer_dict[kmer] += 1
    return kmer_dict


def build_graph(kmer_dict):
    """ Builds Networkx DiGraph according to a kmer dictionary.
    """
    digraph = nx.DiGraph()
    digraph.add_weighted_edges_from(
        [(key[:-1], key[1:], value) for key, value in kmer_dict.items()]
    )
    return digraph


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """ Removes paths.
    """
    for path in path_list:
        if len(path) == 2 and not delete_entry_node and not delete_sink_node:
            graph.remove_edge(path[0], path[1])
        else:
            for i, node in enumerate(path):
                if i == 0 and delete_entry_node:
                    graph.remove_node(node)
                elif i == len(path) - 1 and delete_sink_node:
                    graph.remove_node(node)
                elif i not in [0, len(path) - 1]:
                    preds = [pred for pred in graph.predecessors(node)]
                    if len(preds) > 1:
                        graph.remove_edge(path[i - 1], node)
                        break
                    graph.remove_node(node)
    return graph


def std(data):
    """ Computes standard deviation.
    """
    return statistics.stdev(data)


def select_best_path(
    graph,
    path_list,
    path_length,
    weight_avg_list,
    delete_entry_node=False,
    delete_sink_node=False,
):
    """ Selects best path (highest weight and highest length) and
    removes the others.
    """
    best_path = path_list[0]
    best_path_lenght = path_length[0]
    best_path_weight = weight_avg_list[0]
    paths_to_remove = []

    for i in range(1, len(path_list)):
        if weight_avg_list[i] > best_path_weight:
            paths_to_remove.append(best_path)
            best_path = path_list[i]
            best_path_lenght = path_length[i]
            best_path_weight = weight_avg_list[i]
        elif weight_avg_list[i] == best_path_weight:
            if path_length[i] > best_path_lenght:
                paths_to_remove.append(best_path)
                best_path = path_list[i]
                best_path_lenght = path_length[i]
                best_path_weight = weight_avg_list[i]
            else:
                paths_to_remove.append(path_list[i])
        else:
            paths_to_remove.append(path_list[i])
    return remove_paths(graph, paths_to_remove, delete_entry_node, delete_sink_node)


def path_average_weight(graph, path):
    """ Compute average weight of a given path in a graph.
    """
    return sum(
        [graph[path[i]][path[i + 1]]["weight"] for i in range(len(path) - 1)]
    ) / (len(path) - 1)


def solve_bubble(graph, ancestor_node, descendant_node):
    """ Cleans graph of bubble between two nodes.
    """
    # get all paths between ancestor node and descendant node
    paths = [
        path
        for path in nx.all_simple_paths(
            graph, source=ancestor_node, target=descendant_node
        )
    ]
    if len(paths) > 1:
        weights = [path_average_weight(graph, path) for path in paths]
        lenghts = [len(path) - 1 for path in paths]
        # clean bubble
        graph = select_best_path(graph, paths, lenghts, weights)
    return graph


def search_for_ancestors(graph, currents, nodes):
    """ FIXME: Deprecated. search_for_common_ancestors should be able to replace it.
    """
    preds = []
    for current in currents:
        preds += [pred for pred in graph.predecessors(current)]
    nodes += preds
    if len(preds) < 1:
        return []
    if len(set(nodes)) != len(nodes):
        uniq = []
        duplicate = []
        for node in nodes:
            if node not in uniq:
                uniq.append(node)
            elif node not in duplicate:
                duplicate.append(node)
        return duplicate
    return search_for_ancestors(graph, preds, nodes)


def simplify_bubbles(graph):
    """ Removes all bubbles of a graph 
    """
    nodes_to_check = [node for node, preds in graph.pred.items() if preds]
    nodes_to_check = [
        node
        for node in nodes_to_check
        if len([pred for pred in graph.predecessors(node)]) > 1
    ]

    for node in nodes_to_check:
        ancestors = search_for_ancestors(graph, [node], [node])
        while len(ancestors) > 0:
            anc_des_list = [[anc, node] for anc in ancestors]
            for ancestor, descendant in anc_des_list:
                graph = solve_bubble(graph, ancestor, descendant)
            if len([pred for pred in graph.predecessors(node)]) < 2:
                ancestors = []
            else:
                ancestors = search_for_ancestors(graph, [node], [node])

    return graph


def solve_entry_tips(graph, starting_nodes):
    """ Removes wrong entry paths from a graph
    """
    if len(starting_nodes) > 1:
        ending_nodes = get_sink_nodes(graph)
        for e_node in ending_nodes:
            path_list = [
                [
                    path
                    for path in nx.all_simple_paths(graph, source=s_node, target=e_node)
                ]
                for s_node in starting_nodes
            ]
            path_list = [path[0] for path in path_list if len(path) > 0]
            junctions = [node for node, succs in graph.pred.items() if succs]
            junctions = [
                node
                for node in junctions
                if len([succ for succ in graph.predecessors(node)]) > 1
            ]
            if len(junctions) > 0:
                common_junction = None
                for junction in junctions:
                    check = True
                    index = 0
                    while check and index < len(path_list):
                        if junction not in path_list[index]:
                            check = False
                        index += 1
                    if check:
                        common_junction = junction
                path_list = [
                    path[: path.index(common_junction) + 1] for path in path_list
                ]
                lenghts = [len(path) - 1 for path in path_list]
                weights = [path_average_weight(graph, path) for path in path_list]
                graph = select_best_path(
                    graph, path_list, lenghts, weights, delete_entry_node=True
                )
    return graph


def solve_out_tips(graph, ending_nodes):
    """ Removes wrong output paths from a graph
    """
    if len(ending_nodes) > 1:
        starting_nodes = get_starting_nodes(graph)
        for s_node in starting_nodes:
            path_list = [
                [
                    path
                    for path in nx.all_simple_paths(graph, source=s_node, target=e_node)
                ]
                for e_node in ending_nodes
            ]
            path_list = [path[0] for path in path_list if len(path) > 0]
            junctions = [node for node, succs in graph.succ.items() if succs]
            junctions = [
                node
                for node in junctions
                if len([succ for succ in graph.successors(node)]) > 1
            ]
            if len(junctions) > 0:
                common_junction = None
                for junction in junctions:
                    check = True
                    index = 0
                    while check and index < len(path_list):
                        if junction not in path_list[index]:
                            check = False
                        index += 1
                    if check:
                        common_junction = junction
                path_list = [path[path.index(common_junction) :] for path in path_list]
                lenghts = [len(path) - 1 for path in path_list]
                weights = [path_average_weight(graph, path) for path in path_list]
                graph = select_best_path(
                    graph, path_list, lenghts, weights, delete_sink_node=True
                )

    return graph


def get_starting_nodes(graph):
    """ Returns starting nodes as list.
    """
    return [node for node, preds in graph.pred.items() if not preds]


def get_sink_nodes(graph):
    """ Returns ending nodes as list.
    """
    return [node for node, succs in graph.succ.items() if not succs]


def get_contigs(graph, starting_nodes, ending_nodes):
    """ Returns contigs according to starting and ending nodes.
    """
    contigs = []
    for s_node in starting_nodes:
        for e_node in ending_nodes:
            for path in nx.all_simple_paths(graph, source=s_node, target=e_node):
                contig = "".join([s_node[0]] + [n[0] for n in path[1:-1]] + [e_node])
                contigs.append((contig, len(contig)))
    return contigs


def save_contigs(contigs_list, output_file):
    """ Saves contigs into file.
    """
    def fill(text, width=80):
        """Split text with a line return to respect fasta format"""
        return os.linesep.join(text[i : i + width] for i in range(0, len(text), width))

    message = ""
    for index, contig in enumerate(contigs_list):
        message += fill(f">contig_{index} len={contig[1]}\n")
        message += fill(contig[0] + "\n")

    ofile = open(output_file, "w")
    ofile.write(message)
    ofile.close()


# def draw_graph(graph, graphimg_file):
#     """Draw the graph
#     """
#     fig, ax = plt.subplots()
#     elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
#     #print(elarge)
#     esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
#     #print(elarge)
#     # Draw the graph with networkx
#     #pos=nx.spring_layout(graph)
#     pos = nx.random_layout(graph)
#     nx.draw_networkx_nodes(graph, pos, node_size=6)
#     nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
#     nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5,
#                            edge_color='b', style='dashed')
#     #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
#     # save image
#     plt.savefig(graphimg_file)


# ==============================================================
# Main program
# ==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    kmer_dict = build_kmer_dict(args.fastq_file, args.kmer_size)
    graph = build_graph(kmer_dict)
    graph = simplify_bubbles(graph)
    starting_nodes = get_starting_nodes(graph)
    ending_nodes = get_sink_nodes(graph)
    graph = solve_entry_tips(graph, starting_nodes)
    graph = solve_out_tips(graph, ending_nodes)
    starting_nodes = get_starting_nodes(graph)
    ending_nodes = get_sink_nodes(graph)
    contigs_list = get_contigs(graph, starting_nodes, ending_nodes)
    # draw_graph(graph, "graph.png")
    save_contigs(contigs_list, args.output_file)


if __name__ == "__main__":
    main()
