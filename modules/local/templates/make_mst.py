#!/usr/bin/env python


"""Make MST."""


import csv
import logging
import platform
import pydot
import sys
import yaml
from heapq import heappop, heappush
from pathlib import Path


logger = logging.getLogger()


def make_graph(lo_concat_pairwise_diffs):
    """
    Turns a list of [(G1, G2, V1),(G2, G1, V1),...] tuples into a graph, which
      is in this case a dictionary of dictionaries. G1 and G2 are the Genomes
      (isolates), V1 (Value) is the number of either SNPs or mutation events.
    param: list lo_concat_pairwise_diffs = list of (G1, G2, V1) and its
           inverse, (G2, G1, V1), which are both needed by prim_mst()
    return: a graph (a dict of dict), input for prim_mst()
    """

    graph = {}

    # unpack each tuple in the list
    for pair in lo_concat_pairwise_diffs:
        G1, G2, V1 = pair

        # if first node already in dictionary, add second node and edge
        if G1 in graph.keys():
            graph[G1][G2] = V1
        # new node, add to dictionary
        else:
            graph[G1] = {G2:V1}

    return graph


def prim_mst(graph):
    """
    Construct the Minimum Spanning Tree for a graph and starting node, using
      Prim's greedy algorithm.
    param: dict graph = a dict of dict for nodes and edge lengths
    param: str start = name of the isolate to start the MST search
    return: dict MST = minimum spanning tree for the graph, in the form:
      G1:[G2, G3]
    """

    start = next(iter(graph))
    lo_nodes = []  # nodes in MST
    MST = {}       # the MST

    # Priority Queue (weight , previous_node , current_node)
    queue = [(0, None, start)]

    # queue empties due to heappop command
    while queue:
        # choose edge with smallest weight, which is done automatically
        # by heappop, the 'weight' varianble is never used
        weight, previous_node, current_node = heappop(queue)

        #skip any vertices already in the MST
        if current_node in lo_nodes:
            continue

        # add current node to list
        lo_nodes.append(current_node)

        # add to the MST structure, which is a dictionary node:list-of-nodes
        if previous_node is None:
            pass
        # at leaast one key:value pair is already present, add the new value
        #  to the list
        elif previous_node in MST:
            MST[previous_node].append(current_node)
        # add new key:value pair to MST dictionary
        else:
            MST[previous_node]=[current_node]

        # retrieve new node and edge weight from graph and add to queue
        for new_node, edge_weight in graph[current_node].items():
            heappush(queue, (edge_weight, current_node, new_node))

    return MST


def weighted_mst(lo_concat_pairwise_diffs, MST):
    """
    Add the weights to the MST. The MST is a dict of key:list items, which is
      less ideal for plotting and lacks the number of mutation events or SNPs
      (aka INT). This function converts the dict into a list of tuples and
      adds back the number of mutation events or SNPs.
    param: list lo_concat_pairwise_diffs = [(G1, G2, INT), (G2, G1, INT), ...]
           where G1 and G2 are isolate names and INT the number of mutation
           events or SNPs
    param: dict MST = dictionary of the form G1:[G2, G3] for the MST
    return: list lo_weighted_MST = a list of tuples [(G1, G2, INT),...] to
                                   draw the MST
    """

    lo_MST = []
    lo_weighted_MST = []

    # turn the dict G1:[G2,G3] into a list [(G1,G2), (G1:G3),...]
    for key in MST.keys():
        lo_vals = MST.get(key)
        for val in lo_vals:
            lo_MST.append((key, val))

    # use the list to select those triplet tuples (G1, G2, (V1))
    # that belong to the MST
    for nodes in lo_MST:
        for tupel in lo_concat_pairwise_diffs:
            if nodes[0] == tupel[0] and nodes[1] == tupel[1]:
                lo_weighted_MST.append(tupel)

    # returns an MST with the weight (V1) added
    return lo_weighted_MST


def get_ref_colors(genome, reference):
    """
    Returns a color for the MST if the reference's name is in the list,
      else white.
    param: str ref = name of a reference strain
    """

    color = 'white'

    if genome == 'Lpn':
        do_ref_col = {
            'ATCC_43290':'aliceblue',
            'Corby':'bisque',
            'D-7158':'aquamarine',
            'Dallas_1E':'azure',
            'Detroit-1':'lightblue',
            'F-4185':'khaki',
            'F4468':'yellowgreen',
            'HL06041035':'orange',
            'Lens':'palegreen',
            'Lorraine':'plum',
            'NCTC11286':'lightgrey',
            'Paris':'navajowhite',
            'Philadelphia_1':'lavender',
            'Pontiac':'coral',
            'ST23':'gold',
            'ST42':'dodgerblue',
            'Toronto-2005':'mintcream',
            'pLPP':'beige'}

        color = do_ref_col.get(reference, 'white')

    return color


def draw_graph(reference, color):
    """
    Returns a graph object containing the reference strain to start with.
    param: str ref = name of a reference strain
    param: str color = name of a color for all nodes (except the reference)
    return: graph_object = a graph object with the reference strain as first
            node; e.g.: draw_graph('a_name', 'red') returns:
            graph G {layout=dot; a_name [shape=box, style=filled,
                     fillcolor=red]; }

    """

    # make a new, undirected graph
    # layout 'dot' gives a nice, hirachical layout, alternatives include
    #  'neato', 'fdp', 'twopi', stay away from 'circo'
    graph_object = pydot.Dot(graph_type='graph', layout='dot')

    # add the reference strain
    graph_object.add_node(pydot.Node(reference, shape='box', style='filled',
                                     fillcolor=color))
    logger.info('## draw_graph() complete')
    return graph_object


def adding_nodes(edge, graph_object, color):
    """
    Adds a new node (= genome) to the graph_object and connects it with an
      edge to an existing node.
    param: tuple edge = (G1, G2, V1) for names of genomes 1 and 2 and count of
           mutation events or SNPs
    param: graph_object = the (growing) graph object
    param: str color = name of a color for all nodes (except the reference)
    return: a graph object with added nodes and edges
    """

    # unpack the tuple of (G1, G2, V1)
    G1, G2, V1 = edge

    # add new nodes for G1 and G2 to the graph_object
    graph_object.add_node(pydot.Node(G1, style='filled', fillcolor=color))
    graph_object.add_node(pydot.Node(G2, style='filled', fillcolor=color))

    # connect G1 to G2, add the number of mutation events or SNPs as edge label
    graph_object.add_edge(pydot.Edge(G1, G2, label=V1, arrowhead='none',
                                     color='blue', fontcolor='red'))

    return graph_object


def finish_draw_graph(graph_object, reference, color):
    """
    Changes the colors and shapes of the reference and query in the MST.
    param: graph_object = the graph object
    param: str ref = name of a reference strain
    param: str isolate = isolate name, e.g.: 'IDR001234'
    param: str color = name of a color for all nodes (except the reference)
    return: a graph object with highlighted reference and query
    """

    graph_object.add_node(pydot.Node(reference, shape='box', style='filled',
                                     fillcolor='orange'))

    # graph_object.add_node(pydot.Node(isolate, shape='box', style='filled',
    #                                  fillcolor=color))

    return graph_object


def make_mst(concat_pairwise_diffs_file, mst_file, report_file, genome, suffix, reference):
    """
    main function
    param: str isolate = name of the bacterial isolate, user supplied
    param: str ref_fa_file = name of a reference strain's FASTA file
    param: list lo_concat_pairwise_diffs = differences (mutation events or
           SNPs) between isolate pairs in the same cluster
    param: str FILE_NAME = file name, 'MST_ME.png' or 'MST_SNP.png'
    param: str suffix = 'ME' or 'SNP'
    """

    lo_concat_pairwise_diffs = []

    with open(concat_pairwise_diffs_file, newline='') as infile:
        reader = csv.reader(infile)
        for G1, G2, V1 in reader:
            lo_concat_pairwise_diffs.append((G1, G2, int(V1)))

    # formats the list of data into a graph dict
    graph = make_graph(lo_concat_pairwise_diffs)
    logger.info('## make_graph() completed')

    # returns a Minimum Spanning Tree
    MST = prim_mst(graph)
    logger.info('## prim_mst() completed')

    # converts the MST dict back into a list of tuples
    lo_weighted_MST = weighted_mst(lo_concat_pairwise_diffs, MST)
    logger.info('## write_to_log() completed')

    # sets the background color for the nodes in the graph drawing, will be
    # 'white' if reference is not in the reference:color dict
    color = get_ref_colors(genome, reference)

    # generate a graph_object for drawing, start with the reference strain as
    #  first node
    graph_drawing = draw_graph(reference, color)
    logger.info('## draw_graph() completed')

    # adding nodes and egdes to the graph_object
    for edge in lo_weighted_MST:
        adding_nodes(edge, graph_drawing, color)
    logger.info('## add_node() completed')

    # highlights the reference
    finish_draw_graph(graph_drawing, reference, color)

    # save the drawn graph_object to file
    graph_drawing.write_png(mst_file)
    logger.info('## write_png() completed for', suffix)

    # write note to report file
    with open(report_file, 'a') as report_file:
        print('\\nFigure: Minimum Spanning tree (' + suffix + ')\\n', file=report_file)

    logger.info('## Added a Minimum Spanning Tree (' + suffix + ').')


if __name__ == "__main__":
    logging.basicConfig(filename="$log_file", level="$log_level", format="[%(levelname)s] %(message)s")

    versions = {}
    versions["${task.process}"] = {
        "python": platform.python_version(),
        "yaml": yaml.__version__,
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f, default_flow_style=False)

    sys.exit(make_mst("$concat_pairwise_diffs", "$mst", "$report", "$genome", "$suffix", "$meta.ref"))
