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
import networkx as nx
import matplotlib
from operator import itemgetter
import random
from random import randint
import statistics
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
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()


def read_fastq(fastq_file):
    with open(fastq_file, 'r') as fastq_file:
        for _ in fastq_file:
            yield next(fastq_file)[:-1]  # remove '\n' character
            next(fastq_file)  # skip useless line
            next(fastq_file)  # skip useless line


def cut_kmer(read, kmer_size):
    for kmer in [read[i:i+kmer_size] for i in range(0, len(read)-kmer_size+1)]:
        yield kmer


def build_kmer_dict(fastq_file, kmer_size):
    kmer_dict = {}
    for sequence in read_fastq(fastq_file):
        for kmer in cut_kmer(sequence, kmer_size):
            if kmer not in kmer_dict.keys():
                kmer_dict[kmer] = 1
            else:
                kmer_dict[kmer] += 1
    return kmer_dict

def build_graph(kmer_dict):
    digraph = nx.DiGraph()
    digraph.add_weighted_edges_from([(key[:-1], key[1:], value) for key, value in kmer_dict.items()])
    return digraph


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    pass

def std(data):
    pass


def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    pass

def path_average_weight(graph, path):
    pass

def solve_bubble(graph, ancestor_node, descendant_node):
    pass

def simplify_bubbles(graph):
    pass

def solve_entry_tips(graph, starting_nodes):
    pass

def solve_out_tips(graph, ending_nodes):
    pass

def get_starting_nodes(graph):
    return [node for node, preds in graph.pred.items() if not preds]

def get_sink_nodes(graph):
    return [node for node, succs in graph.succ.items() if not succs]

def get_contigs(graph, starting_nodes, ending_nodes):
    contigs = []
    for s_node in starting_nodes:
        for e_node in ending_nodes:
            for path in nx.all_simple_paths(graph, source=s_node, target=e_node):
                contig = "".join([s_node[0]] + [n[0] for n in path[1:-1]] + [e_node])
                contigs.append((contig, len(contig)))
    return contigs

def save_contigs(contigs_list, output_file):
    def fill(text, width=80):
        """Split text with a line return to respect fasta format"""
        return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))
    message = ""
    for index, contig in enumerate(contigs_list):
        message += fill(f">contig_{index} len={contig[1]}\n")
        message += fill(contig[0]+"\n")

    ofile = open(output_file, "w")
    ofile.write(message)
    ofile.close()


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    kmer_dict = build_kmer_dict(args.fastq_file, args.kmer_size)
    graph = build_graph(kmer_dict)
    starting_nodes = get_starting_nodes(graph)
    ending_nodes = get_sink_nodes(graph)
    contigs_list = get_contigs(graph, starting_nodes, ending_nodes)
    save_contigs(contigs_list, args.output_file)


if __name__ == '__main__':
    main()
