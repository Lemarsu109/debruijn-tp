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
#import matplotlib as plt
from operator import itemgetter
import random
from networkx.algorithms.shortest_paths.unweighted import predecessor

from networkx.readwrite import gpickle
random.seed(9001)
from random import randint
import statistics

__author__ = "Aidan Bonnefond"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Aidan Bonnefond"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Aidan Bonnefond"
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
                        default=22, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    with open(fastq_file, "rt") as f:
        for i,line in enumerate(f):
            #print(line)
            if i % 4 == 1:
                #print(line)
                yield line


def cut_kmer(read, kmer_size):
    for i in range(len(read)-kmer_size):
        #print(read[i:i+kmer_size])
        yield read[i:i+kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
    kmer_dict = {}
    gen_line = read_fastq(fastq_file)
    for line in gen_line:
        gen_kmer = cut_kmer(line, kmer_size)
        for kmer in gen_kmer:
            if kmer in kmer_dict:
                kmer_dict[kmer]+=1
            else:
                kmer_dict[kmer]=1
    return kmer_dict


def build_graph(kmer_dict):
    graph = nx.DiGraph()
    for i in kmer_dict:
        #print("prefix : "+str(i[0:len(i)-1])+" suffix : "+str(i[1:len(i)])+" weigth : "+str(kmer_dict[i]))
        graph.add_edge(i[:-1], i[1:], weight=kmer_dict[i])
    return graph


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    for path in path_list:
        if delete_entry_node and delete_sink_node:
            graph.remove_nodes_from(path)
        elif delete_entry_node == False and delete_sink_node:
            graph.remove_nodes_from(path[1:])
        elif delete_entry_node and delete_sink_node == False:
            graph.remove_nodes_from(path[:-1])
        else:
            graph.remove_nodes_from(path[1:-1])

def std(data):
    return statistics.stdev(data)


def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    if std(weight_avg_list) > 0:
        for i in path_list:
            if weight_avg_list[i] != max(weight_avg_list):
                remove_paths(graph, path_list[i], delete_entry_node, delete_sink_node)
        return graph
    elif std(path_length) > 0:
        for i in path_list:
            if path_length[i] != max(path_length):
                remove_paths(graph, path_list[i], delete_entry_node, delete_sink_node)
        return graph
    else:
        rand = randint(0, len(path_list))
        for i in path_list:
            if i != rand:
                    remove_paths(graph, path_list[i], delete_entry_node, delete_sink_node)
        return graph


def path_average_weight(graph, path):
    av_mean = 0
    for u,v,d in graph.subgraph(path).edges(data=True):
        av_mean = statistics.mean([d["weight"] for u,v,d in graph.subgraph(path).edges(data=True)])
    return av_mean

def solve_bubble(graph, ancestor_node, descendant_node):
    path_list = nx.all_simple_paths(graph, ancestor_node, descendant_node)
    path_length = []
    weight_avg_list = []
    for i in path_list:
        path_length.append(len(path_list[i]))
        weight_avg_list.append(path_average_weight(graph, path_list[i]))
    select_best_path(graph, path_list, path_length, weight_avg_list)

def simplify_bubbles(graph):
    bubble = False
    for n in graph.nodes():
            list_pred = graph.predecessors(n)
            if len(list(list_pred)) > 1:
                for i in graph.predecessors(n):
                    for j in graph.predecessors(n):
                        if i !=j:
                            ancetre = nx.lowest_common_ancestor(graph, i, j)
                            if ancetre != None:
                                bubble = True
                    if bubble:
                        break
                if bubble:
                    break
            if bubble:
                break
    if bubble:
        graph = simplify_bubbles(solve_bubble(graph, ancetre, n))
    return graph

def solve_entry_tips(graph, starting_nodes):
    path_list = []
    weight_avg_list = []
    for n in graph.nodes():
        i = 0
        for p in graph.predecessors(n):
            P = []
            P.append(p)
            path_list.append(P)
            weight_avg_list.append(path_average_weight(graph, p))
            for s in starting_nodes:
                if p == s:
                    i += 1
        if i > 1:
            graph = solve_entry_tips(select_best_path(graph, path_list, weight_avg_list))
    return graph
    

def solve_out_tips(graph, ending_nodes):
    path_list = []
    weight_avg_list = []
    for n in graph.nodes():
        i = 0
        for p in graph.successors(n):
            P = []
            P.append(p)
            path_list.append(P)
            weight_avg_list.append(path_average_weight(graph, p))
            for s in ending_nodes:
                if p == s:
                    i += 1
        if i > 1:
            graph = solve_entry_tips(select_best_path(graph, path_list, weight_avg_list))
    return graph

def get_starting_nodes(graph):
    for n in graph.nodes():
        pred = [p for p in graph.predecessors(n)]
    return pred

def get_sink_nodes(graph):
    for n in graph.nodes():
        succ = [p for p in graph.successors(n)]
    return succ

def get_contigs(graph, starting_nodes, ending_nodes):
    contig = []
    for i in starting_nodes:
        for j in ending_nodes:
            for path in nx.all_simple_paths(graph, i, j):
                seq = path[0]
                for s in path:
                    seq+=s[-1]
                contig.append((seq, len(seq)))
    return contig
                

def save_contigs(contigs_list, output_file):
    with open(output_file, "w") as o:
        for i,contig in enumerate(contigs_list):
            o.write(">contig_"+str(i)+" len="+str(contig[1])+"\n")
            o.write(fill(contig[0])+"\n")
    o.close()


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def draw_graph(graph, graphimg_file):
    """Draw the graph
    """                                    
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


def save_graph(graph, graph_file):
    """Save the graph with pickle
    """
    with open(graph_file, "wt") as save:
            gpickle.dump(graph, save)


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    k = 22
    graph = build_graph(build_kmer_dict(args.fastq_file, k))
    #print(graph)
    print(get_starting_nodes(graph),get_sink_nodes(graph))
    graph = simplify_bubbles(graph)
    #print(graph)
    graph = solve_entry_tips(graph, get_starting_nodes(graph))
    graph = solve_out_tips(graph, get_sink_nodes(graph))
    #print(graph)
    contig = get_contigs(graph, get_starting_nodes(graph),get_sink_nodes(graph))
    #print(contig)
    #print(get_starting_nodes(graph),get_sink_nodes(graph))
    save_contigs(contig, "/Users/aidanbonnefond/debruijn-tp/data/eva71_hundred_reads.fq")


    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)
    # Save the graph in file
    # if args.graph_file:
    #     save_graph(graph, args.graph_file)


if __name__ == '__main__':
    main()
