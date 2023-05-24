#   #!/usr/bin/env python
#   -*- coding: utf-8 -*-
#   ******************************************************************************
#     Copyright (c) 2023.
#     Developed by Yifei Lu
#     Last change on 4/5/23, 9:19 AM
#     Last change by yifei
#    *****************************************************************************
from collections import Counter
import networkx as nx
from ..network import Network


def create_graph(network: Network):
    pipelines = network.pipelines
    l_edges = []
    for pipeline in pipelines.values():
        l_edges.append({pipeline.inlet_index, pipeline.outlet_index})
    G = nx.Graph()
    G.add_edges_from(l_edges)
    return G


def graph_nodal_degree_counter(G: nx.Graph):
    degrees = [d for (i, d) in G.degree()]
    return dict(sorted(Counter(degrees).items()))


def nodes_with_degree_n(G: nx.Graph, n: int):
    degree_counter = graph_nodal_degree_counter(G)
    if n in degree_counter.keys():
        return degree_counter[n]
    else:
        print(f"There is no node in this network with degree {n}!")
        return None
