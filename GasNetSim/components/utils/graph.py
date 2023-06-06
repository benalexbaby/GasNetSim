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
    connections = network.connections
    l_edges = []
    for connection in connections.values():
        l_edges.append({connection.inlet_index, connection.outlet_index})
    G = nx.Graph()
    G.add_edges_from(l_edges)
    return G


def graph_nodal_degree_counter(G: nx.Graph):
    degrees = [d for (i, d) in G.degree()]
    return dict(sorted(Counter(degrees).items()))


def nodes_with_degree_n(graph, n):
    node_indices = []
    for node in graph.nodes():
        degree = graph.degree(node)
        if degree == n:
            node_indices.append(node)
    return node_indices
