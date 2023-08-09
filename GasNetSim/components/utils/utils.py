#   #!/usr/bin/env python
#   -*- coding: utf-8 -*-
#   ******************************************************************************
#     Copyright (c) 2022.
#     Developed by Yifei Lu
#     Last change on 3/28/22, 12:44 AM
#     Last change by yifei
#    *****************************************************************************
import math
from pyparsing import col
from collections import OrderedDict
from scipy import sparse
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def create_connection_matrix(n_nodes: int, components: dict, component_type: int, sparse_matrix: bool = False):
    row_ind = list()
    col_ind = list()
    data = list()

    if not sparse_matrix:
        cnx = np.zeros((n_nodes, n_nodes))

    for comp in components.values():
        i = comp.inlet_index - 1
        j = comp.outlet_index - 1
        if sparse_matrix:
            row_ind.append(i)
            col_ind.append(j)
            data.append(component_type)
        else:
            cnx[i][j] = component_type
            cnx[j][i] = component_type

    if sparse_matrix:
        cnx = sparse.coo_matrix((data, (row_ind, col_ind)))
    return cnx


def levenberg_marquardt_damping_factor(m, s, b):
    return 10 ** (m * math.log10(s + b))


def delete_matrix_rows_and_columns(matrix, to_remove):
    new_matrix = matrix

    new_matrix = np.delete(new_matrix, to_remove, 0)  # delete rows
    new_matrix = np.delete(new_matrix, to_remove, 1)  # delete columns

    return new_matrix


def jacobian_matrix_condition_number(matrix):
    print(f"The condition number of the matrix is {np.linalg.cond(matrix)}.")


def print_n_largest_absolute_values(n, values):
    sorted_values = sorted([abs(x) for x in values])
    print(sorted_values[-n::-1])
    return None


def calculate_nodal_inflow_states(nodes, connections, mapping_connections, flow_matrix):
    nodal_total_inflow = np.sum(np.where(flow_matrix > 0, flow_matrix, 0), axis=1)

    nodal_gas_inflow_composition = dict()
    nodal_gas_inflow_temperature = dict()

    for i_node, node in nodes.items():  # iterate over all nodes
        inflow_from_node = np.where(flow_matrix[i_node-1] > 0)[0]  # find the supplying nodes
        if len(inflow_from_node) == 0:
            pass
        else:
            inflow_from_node += 1

        total_inflow_comp = dict()
        total_inflow = nodal_total_inflow[i_node-1]
        total_inflow_temperature_times_flow_rate = 0

        for inlet_index in inflow_from_node:
            gas_composition = nodes[inlet_index].gas_mixture.composition
            connections[mapping_connections[i_node - 1][inlet_index - 1]].gas_mixture.composition = gas_composition
            inflow_rate = flow_matrix[i_node-1][inlet_index-1]
            inflow_temperature = connections[mapping_connections[i_node-1][inlet_index-1]].calc_pipe_outlet_temp()

            # Sum up flow rate * temperature
            total_inflow_temperature_times_flow_rate += inflow_rate * inflow_temperature

            # create a OrderedDict to store gas flow fractions
            gas_flow_comp = OrderedDict({gas: comp * inflow_rate for gas, comp in gas_composition.items()})
            for gas, comp in gas_flow_comp.items():
                if total_inflow_comp.get(gas) is None:
                    total_inflow_comp[gas] = comp
                else:
                    total_inflow_comp[gas] += comp

        nodal_gas_inflow_composition[i_node] = {k: v / total_inflow for k, v in total_inflow_comp.items()}

        if total_inflow != .0:
            nodal_gas_inflow_temperature[i_node] = total_inflow_temperature_times_flow_rate / total_inflow
        else:
            nodal_gas_inflow_temperature[i_node] = np.nan

    return nodal_gas_inflow_composition, nodal_gas_inflow_temperature


def calculate_flow_matrix(network, pressure_bar):
    connections = network.connections
    nodes = network.nodes
    n_nodes = len(nodes)
    flow_mat = np.zeros((n_nodes, n_nodes), dtype=float)

    pressure_index = 0
    for node in nodes.values():
        if node.index not in network.non_junction_nodes:
            node.pressure = pressure_bar[pressure_index] * 1e5
            pressure_index += 1

    for connection in connections.values():
        i = connection.inlet_index - 1
        j = connection.outlet_index - 1
        connection.inlet = nodes[i+1]
        connection.outlet = nodes[j+1]

        flow_direction = connection.determine_flow_direction()

        p1 = nodes[i+1].pressure
        p2 = nodes[j+1].pressure

        slope_correction = connection.calc_pipe_slope_correction()
        temp = connection.calculate_coefficient_for_iteration()

        flow_rate = flow_direction * abs(p1 ** 2 - p2 ** 2 - slope_correction) ** (1 / 2) * temp

        flow_mat[i][j] = - flow_rate
        flow_mat[j][i] = flow_rate

    return flow_mat


def calculate_flow_vector(network, pressure_bar, target_flow):
    flow_matrix = calculate_flow_matrix(network, pressure_bar)
    n_nodes = len(network.nodes.values())
    nodal_flow = np.dot(flow_matrix, np.ones(n_nodes))
    nodal_flow = [nodal_flow[i] for i in range(len(nodal_flow)) if i + 1 not in network.non_junction_nodes]
    delta_flow = target_flow - nodal_flow

    # delta_flow = [delta_flow[i] for i in range(len(delta_flow)) if i + 1 not in network.non_junction_nodes]
    return delta_flow


def plot_network_demand_distribution(network):
    nodes = network.nodes.values()
    node_demand = [n.volumetric_flow for n in nodes if n.volumetric_flow is not None]
    sns.histplot(data=node_demand, stat="probability")
    plt.xlim((min(node_demand)-10, max(node_demand) + 10))
    plt.xlabel("Nodal volumetric flow demand [sm^3/s]")
    plt.show()
    return None


def check_square_matrix(a):
    return a.shape[0] == a.shape[1]

def check_symmetric(a, rtol=1e-05, atol=1e-08):
    return np.allclose(a, a.T, rtol=rtol, atol=atol)


def check_all_off_diagonal_elements(a, criterion):
    res = True

    if check_square_matrix(a):
        pass
    else:
        print("Matrix is not a square matrix!")

    for i in range(a.shape[0]):
        for j in range(a.shape[1]):
            if i != j:
                if criterion == "zero":
                    res = (a[i][j] == 0)
                elif criterion == "positive":
                    res = (a[i][j] > 0)
                elif criterion == "non-negative":
                    res = (a[i][j] >= 0)
                elif criterion == "negative":
                    res = (a[i][j] < 0)
                elif criterion == "non-positive":
                    res = (a[i][j] <= 0)
                else:
                    print("Check the given criterion!")
                    return False
                if res == False:
                    return False
    return res


def check_all_diagonal_elements(a, criterion):
    res = True

    if check_square_matrix(a):
        pass
    else:
        print("Matrix is not a square matrix!")

    if criterion == "zero":
        res = (np.diagonal(a) == 0).all()
    elif criterion == "positive":
        res = (np.diagonal(a) > 0).all()
    elif criterion == "non-negative":
        res = (np.diagonal(a) >= 0).all()
    elif criterion == "negative":
        res = (np.diagonal(a) < 0).all()
    elif criterion == "non-positive":
        res = (np.diagonal(a) <= 0).all()
    else:
        print("Check the given criterion!")
        return False

    return res
