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


def graph_info(G: nx.Graph):
    degrees = []
    for (i, d) in G.degree():
        degrees.append(d)
    print(sorted(Counter(degrees).items()))
