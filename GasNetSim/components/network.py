#  #!/usr/bin/env python
#  -*- coding: utf-8 -*-
#  ******************************************************************************
#    Copyright (c) 2022.
#    Developed by Yifei Lu
#    Last change on 1/17/22, 11:21 AM
#    Last change by yifei
#   *****************************************************************************

import numpy as np
from typing import Tuple
from scipy import sparse
import logging
import copy
from collections import OrderedDict

from .utils.gas_mixture.heating_value import *
from .utils.utils import *
from .node import *
from .pipeline import *
from .utils import *

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class Network:
    """
    Network class
    """

    def __init__(self, nodes: dict, pipelines=None, compressors=None, resistances=None, shortpipes=None):
        """

        :param nodes:
        :param pipelines:
        """
        self.nodes = nodes
        self.pipelines = pipelines
        self.compressors = compressors
        self.resistances = resistances
        self.shortpipes = shortpipes
        self.connections = self.all_edge_components()
        self.connection_matrix = self.create_connection_matrix()
        self.reference_nodes = self.find_reference_nodes()
        self.non_junction_nodes = self.find_non_junction_nodes()
        self.junction_nodes = self.find_junction_nodes()

    def all_edge_components(self):
        connections = dict()

        all_edge_classes = [self.pipelines, self.resistances, self.compressors, self.shortpipes]

        i_connection = 0

        for edge_class in all_edge_classes:
            if edge_class is not None:
                for edge in edge_class.values():
                    connections[i_connection] = edge
                    i_connection += 1

        return connections

    def mapping_of_connections(self):
        n_nodes = len(self.nodes.values())
        mapping = np.zeros((n_nodes, n_nodes))
        for i_connection, connection in self.connections.items():
            i = connection.inlet_index - 1
            j = connection.outlet_index - 1
            mapping[i][j] = mapping[j][i] = i_connection
        return mapping

    def find_reference_nodes(self) -> list:
        """
        Find reference nodes, where pressure are pre-set
        :return: List of Node class of pressure-referenced nodes
        """
        pressure_ref_nodes = list()

        for node in self.nodes.values():
            if node.pressure is not None:
                pressure_ref_nodes.append(node.index)

        return pressure_ref_nodes

    # def find_ref_nodes_index(self):
    #     """
    #     To get indices of the referenced nodes
    #     :return: List of pressure-referenced nodes indices and list of temperature-referenced nodes indices
    #     """
    #     pressure_ref_nodes_index = list()
    #     temperature_ref_nodes_index = list()
    #
    #     for node_index, node in self.nodes.items():
    #         if node.pressure is not None:
    #             pressure_ref_nodes_index.append(node_index)
    #         if node.temperature is not None:
    #             temperature_ref_nodes_index.append(node_index)
    #
    #     return pressure_ref_nodes_index, temperature_ref_nodes_index

    def demand_nodes_supply_pipelines(self):
        """
        Find supply pipelines for the demand nodes
        :return:
        """
        nodal_supply_pipelines = dict()

        if self.connections is not None:
            for i_connection, connection in self.connections.items():
                if connection.flow_rate is None or connection.flow_rate > 0:
                    if nodal_supply_pipelines.get(connection.outlet_index) is not None:
                        nodal_supply_pipelines[connection.outlet_index].append(i_connection)
                    else:
                        nodal_supply_pipelines[connection.outlet_index] = [i_connection]
                elif connection.flow_rate < 0:
                    if nodal_supply_pipelines.get(connection.inlet_index) is not None:
                        nodal_supply_pipelines[connection.inlet_index].append(i_connection)
                    else:
                        nodal_supply_pipelines[connection.inlet_index] = [i_connection]

        return OrderedDict(sorted(nodal_supply_pipelines.items()))

    def convert_energy_flow_to_volumetric_flow(self, base='HHV'):
        for node in self.nodes.values():
            gas_comp = node.get_mole_fraction()
            standard_density = GasMixture(pressure=101325, temperature=288.15, composition=gas_comp).density
            LHV, HHV = calc_heating_value(node.gas_mixture)
            if base == 'HHV':
                heating_value = HHV/1e6*standard_density  # MJ/sm3
            elif base == 'LHV':
                heating_value = LHV/1e6*standard_density  # MJ/sm3
            else:
                raise ValueError
            if node.flow is not None:
                node.flow /= heating_value
                # try:
                #     h2_fraction = node.gas_mixture.zs[node.gas_mixture.components.index('hydrogen')]
                # except:
                #     h2_fraction = 0
                # node.flow /= (h2_fraction * 12.09 + (1-h2_fraction) * 38.28)
        return None

    def create_connection_matrix(self, sparse_matrix=False):
        # TODO change the index number
        n_nodes = len(self.nodes.values())
        pipelines = self.pipelines
        compressors = self.compressors
        resistances = self.resistances
        shortpipes = self.shortpipes
        if sparse_matrix:
            row_ind = list()
            col_ind = list()
            data = list()

        # Build a matrix to show the connection between nodes
        connection = np.zeros((n_nodes, n_nodes))

        if pipelines is not None:
            for pipe in pipelines.values():
                i = pipe.inlet_index - 1
                j = pipe.outlet_index - 1
                if sparse_matrix:
                    row_ind.append(i)
                    col_ind.append(j)
                    data.append(1)
                else:
                    connection[i][j] = 1
                    connection[j][i] = 1

        if compressors is not None:
            for compressor in compressors.values():
                i = compressor.inlet_index - 1
                j = compressor.outlet_index - 1
                if sparse_matrix:
                    row_ind.append(i)
                    col_ind.append(j)
                    data.append(2)
                else:
                    connection[i][j] = 2
                    connection[j][i] = 2

        if resistances is not None:
            for resistance in resistances.values():
                i = resistance.inlet_index - 1
                j = resistance.outlet_index - 1
                if sparse_matrix:
                    row_ind.append(i)
                    col_ind.append(j)
                    data.append(3)
                else:
                    connection[i][j] = 3
                    connection[j][i] = 3

        if shortpipes is not None:
            for sp in shortpipes.values():
                i = sp.inlet_index - 1
                j = sp.outlet_index - 1
                if sparse_matrix:
                    row_ind.append(i)
                    col_ind.append(j)
                    data.append(4)
                else:
                    connection[i][j] = 4
                    connection[j][i] = 4

        return connection

    def pressure_initialization(self):
        nodes = self.nodes

        # create a list to store all component resistance
        resistance = list()
        if self.pipelines is not None:
            pipeline_resistance = [[x.inlet_index, x.outlet_index, x.resistance, x.outlet.flow]
                                   for x in self.pipelines.values()]
            resistance += pipeline_resistance
        if self.resistances is not None:
            resistance_resistance = [[x.inlet_index, x.outlet_index, x.resistance, x.outlet.flow]
                                     for x in self.resistances.values()]
            resistance += resistance_resistance
        if self.shortpipes is not None:
            shortpipe_resistance = [[x.inlet_index, x.outlet_index, 0, -x.inlet.flow]
                                    for x in self.shortpipes.values()]
            resistance += shortpipe_resistance

        max_resistance = max([x[2] for x in resistance])
        max_flow = max([x.flow for x in nodes.values() if x.flow is not None])
        pressure_init = [node.pressure for node in nodes.values()]
        # pipeline_with_missing_pressure = copy.deepcopy(pipelines)
        pressure_init_old = list()

        while pressure_init != pressure_init_old:
            pressure_init_old = copy.deepcopy(pressure_init)
            # pipeline_initialized = list()
            for r in resistance:
                i = r[0] - 1  # inlet index
                j = r[1] - 1  # outlet index
                res = r[2]  # resistance
                flow = r[3]
                if pressure_init[i] is None and pressure_init[j] is None:
                    pass
                elif pressure_init[j] is None or pressure_init[i] == pressure_init[j]:
                    pressure_init[j] = pressure_init[i] * (1 - 0.05 * (res / max_resistance) * (flow / max_flow))
                    # if res/max_resistance < 0.001:
                    #     pressure_init[j] = pressure_init[i] * 0.999999
                    # else:
                    #     pressure_init[j] = pressure_init[i] * (1 - 0.05 * (res/max_resistance) * (flow/max_flow))
                        # pressure_init[j] = pressure_init[i] * 0.98
                # elif pressure_init[j] is not None and pressure_init[i] is not None:
                #     if res/max_resistance < 0.001:
                #         pressure_init[j] = min(pressure_init[j], pressure_init[i] * 0.99999)
                #     else:
                #         pressure_init[j] = min(pressure_init[j],
                #                                pressure_init[i] * (1 - 0.05 * (res/max_resistance) * (flow/max_flow)))
                #         # pressure_init[j] = min(pressure_init[j], pressure_init[i] * 0.98)
                elif pressure_init[i] is None and pressure_init[j] is not None:
                    pressure_init[i] = pressure_init[j] / (1 - 0.05 * (res / max_resistance) * (flow / max_flow))
                    # if res/max_resistance < 0.001:
                    #     pressure_init[i] = pressure_init[j] / 0.99999
                    # else:
                    #     pressure_init[i] = pressure_init[j] / (1 - 0.05 * (res/max_resistance) * (flow /max_flow))
                        # pressure_init[i] = pressure_init[j] / 0.98

        return pressure_init

    def newton_raphson_initialization(self):
        """
        Initialization for NR-solver, where the nodal pressures are initialized as 0.98 of inlet pressure and nodal
        temperatures are the same as pipe surrounding temperatures
        :return: Network initial conditions for NR-solver
        """

        nodes = self.nodes
        pipelines = self.pipelines
        resistances = self.resistances

        n_nodes = len(nodes)

        # Build a matrix to show the connection between nodes
        connection = self.create_connection_matrix()

        # TODO consider the case where ref_nodes do not start with index 0
        p_ref_nodes = self.reference_nodes

        for node in self.nodes.values():
            HHV = calc_heating_value(node.gas_mixture)
            if node.flow_type == 'volumetric':
                pass
            elif node.flow_type == 'energy':
                gas_comp = node.get_mole_fraction()
                node.flow = node.flow / HHV * 1e6 / GasMixture(composition=gas_comp,
                                                               temperature=288.15,
                                                               pressure=101325).density
                logging.debug(node.flow)
                node.flow_type = 'volumetric'
            else:
                raise AttributeError(f'Unknown flow type {node.flow_type}!')
        nodal_flow_init = [x.flow if x.flow is not None else 0 for x in nodes.values()]
        pressure_init = [x.pressure for x in nodes.values()]

        # TODO use ambient temperature to initialize outlet temperature
        temperature_init = [x.temperature if x.temperature is not None
                            else 288.15 for x in self.nodes.values()]

        total_flow = sum([x for x in nodal_flow_init if x is not None])

        for n in p_ref_nodes:
            nodal_flow_init[n - 1] = - total_flow / len(p_ref_nodes)

        if None in pressure_init:
            pressure_init = self.pressure_initialization()

        for i in range(len(nodal_flow_init)):
            # TODO change to number of non-reference nodes
            nodes[i + 1].pressure = pressure_init[i]
            nodes[i + 1].flow = nodal_flow_init[i]
            nodes[i + 1].volumetric_flow = nodal_flow_init[i]
            nodes[i + 1].convert_volumetric_to_energy_flow()
            nodes[i + 1].temperature = temperature_init[i]

        if pipelines is not None:
            for index, pipe in pipelines.items():
                pipe.inlet = nodes[pipe.inlet_index]
                pipe.outlet = nodes[pipe.outlet_index]

        if resistances is not None:
            for index, r in resistances.items():
                r.inlet = nodes[r.inlet_index]
                r.outlet = nodes[r.outlet_index]

        return nodal_flow_init, pressure_init, temperature_init

    def jacobian_matrix(self):

        connections = self.connections
        nodes = self.nodes

        reference_nodes = [x-1 for x in self.reference_nodes]  # indices of reference nodes
        non_junction_nodes = [x-1 for x in self.non_junction_nodes]

        n_nodes = len(nodes)

        n_junction_nodes = len(self.junction_nodes)

        jacobian_mat = np.zeros((n_nodes, n_nodes), dtype=np.float)
        flow_mat = np.zeros((n_nodes, n_nodes), dtype=np.float)

        for connection in connections.values():
            i = connection.inlet_index - 1
            j = connection.outlet_index - 1

            flow_mat[i][j] = - connection.calc_flow_rate()
            flow_mat[j][i] = connection.calc_flow_rate()

            if type(connection) is not ShortPipe:
                slope_corr = connection.calc_pipe_slope_correction()
                p1 = connection.inlet.pressure
                p2 = connection.outlet.pressure
                pipeline_coefficient = connection.calculate_coefficient_for_iteration()
                temp_var = (abs(p1 ** 2 - p2 ** 2 - slope_corr)) ** (-0.5)

                if i not in non_junction_nodes and j not in non_junction_nodes:
                    jacobian_mat[i][j] = pipeline_coefficient * p2 * temp_var
                    jacobian_mat[j][i] = pipeline_coefficient * p1 * temp_var
                if i not in non_junction_nodes:
                    jacobian_mat[i][i] += - pipeline_coefficient * p1 * temp_var
                if j not in non_junction_nodes:
                    jacobian_mat[j][j] += - pipeline_coefficient * p2 * temp_var

        jacobian_mat = delete_matrix_rows_and_columns(jacobian_mat, non_junction_nodes)
        # flow_mat = delete_matrix_rows_and_columns(flow_mat, non_junction_nodes)

        return jacobian_mat, flow_mat

    def find_non_junction_nodes(self):
        non_junction_nodes = list()

        if self.shortpipes is not None:
            for sp in self.shortpipes.values():
                non_junction_nodes.append(sp.inlet_index)
        for node in self.nodes.values():
            if node.node_type == 'reference':
                non_junction_nodes.append(node.index)

        return sorted(non_junction_nodes)

    def find_junction_nodes(self):
        return [node.index for node in self.nodes.values() if node.index not in self.non_junction_nodes]

    def newton_raphson_iteration(self, target):
        jacobian_matrix, flow_matrix = self.jacobian_matrix()
        number_of_junction_nodes = len(jacobian_matrix)
        j_mat_inv = np.linalg.inv(jacobian_matrix)

        # Flow balance of connecting pipeline elements at all nodes
        delta_flow = target - np.dot(flow_matrix, np.ones(number_of_junction_nodes))

        # Select only the flows at junction nodes
        delta_flow = [delta_flow[i] for i in range(len(delta_flow)) if i + 1 not in self.junction_nodes]

    def calculate_nodal_inflow_composition(self):
        pass

    def simulation(self):
        logging.debug([x.flow for x in self.nodes.values()])
        # ref_nodes = self.p_ref_nodes_index

        n_nodes = len(self.nodes.keys())
        n_non_junction_nodes = len(self.non_junction_nodes)
        connection_matrix = self.connection_matrix

        init_f, init_p, init_t = self.newton_raphson_initialization()

        max_iter = 100
        n_iter = 0
        # n_non_ref_nodes = n_nodes - len(ref_nodes)

        f = np.array(init_f)
        p = np.array(init_p)
        t = np.array(init_t)
        logging.info(f'Initial pressure: {p}')
        logging.info(f'Initial flow: {f}')

        for i in range(len(init_f)):
            # TODO change to number of non-reference nodes
            self.nodes[i + 1].pressure = init_p[i]
            self.nodes[i + 1].flow = init_f[i]
            self.nodes[i + 1].temperature = init_t[i]
            self.nodes[i + 1].update_gas_mixture()

        if self.pipelines is not None:
            for index, pipe in self.pipelines.items():
                pipe.inlet = self.nodes[pipe.inlet_index]
                pipe.outlet = self.nodes[pipe.outlet_index]
                pipe.update_gas_mixture()

        if self.resistances is not None:
            for index, r in self.resistances.items():
                r.inlet = self.nodes[r.inlet_index]
                r.outlet = self.nodes[r.outlet_index]
                r.update_gas_mixture()

        delta_flow = 0

        record = list()

        while n_iter <= max_iter:
            j_mat, f_mat = self.jacobian_matrix()
            mapping_connections = self.mapping_of_connections()
            nodal_gas_inflow_composition, nodal_gas_inflow_temperature = calculate_nodal_inflow_states(self.nodes, self.connections,
                                          mapping_connections, f_mat)

            for i_node, node in self.nodes.items():
                if nodal_gas_inflow_composition[i_node] == {}:
                    pass
                else:
                    try:
                        self.nodes[i_node].gas_mixture = GasMixture(composition=nodal_gas_inflow_composition[i_node],
                                                                    temperature=self.nodes[i_node].temperature,
                                                                    pressure=self.nodes[i_node].pressure)
                    except Exception:
                        logging.warning(i_node)
                        logging.warning(nodal_gas_inflow_composition[i_node])

            nodal_flow = np.dot(f_mat, np.ones(n_nodes))

            delta_flow = f - nodal_flow

            delta_flow = [delta_flow[i] for i in range(len(delta_flow)) if i + 1 not in self.non_junction_nodes]

            # Update volumetric flow rate target
            for n in self.nodes.values():
                n.convert_energy_to_volumetric_flow()
            f = np.array([x.volumetric_flow if x.flow is not None else 0 for x in self.nodes.values()])

            delta_p = np.linalg.solve(j_mat, delta_flow)  # np.linalg.solve() uses LU decomposition as default
            delta_p /= 2  # divided by 2 to ensure better convergence
            logging.debug(delta_p)

            for i in self.non_junction_nodes:
                delta_p = np.insert(delta_p, i-1, 0)

            p += delta_p  # update nodal pressure list

            # TODO Newton-Raphson with Jacobian adjustment

            for i in self.nodes.keys():
                if i not in self.reference_nodes:
                    self.nodes[i].pressure = p[i - 1]  # update nodal pressure

            for i_connection, connection in self.connections.items():
                connection.inlet = self.nodes[connection.inlet_index]
                connection.outlet = self.nodes[connection.outlet_index]

            record.append(delta_p)

            n_iter += 1
            # j_mat, f_mat = self.jacobian_matrix()
            # nodal_flow = np.dot(f_mat, np.ones(n_nodes))
            # delta_flow = f - nodal_flow
            #
            # delta_flow = np.array([delta_flow[i] for i in range(len(delta_flow)) if i + 1 not in self.non_junction_nodes])
            target_flow = np.array([f[i] for i in range(len(f)) if i + 1 not in self.non_junction_nodes])

            logging.debug(max([abs(x) for x in (delta_flow/target_flow)]))
            logging.debug(delta_p)

            if max([abs(x) for x in (delta_flow/target_flow)]) <= 0.001:
                logger.info(f'Simulation converges in {n_iter} iterations.')
                logger.info(p)
                logger.debug(init_p)
                pipe_h2_fraction = list()

                for i_node in self.non_junction_nodes:
                    self.nodes[i_node].flow = nodal_flow[i_node-1]

                # output connection
                for i_connection, connection in self.connections.items():
                    logger.debug(f'Pipeline index: {i_connection}')
                    logger.debug(f'Pipeline flow rate: {connection.flow_rate}')
                    logger.debug(f'Gas mixture composition: {connection.gas_mixture.composition}')
                    try:
                        pipe_h2_fraction.append(connection.gas_mixture.composition['hydrogen'] * 100)
                    except KeyError:
                        pipe_h2_fraction.append(0)
                logging.debug(pipe_h2_fraction)
                for connection in self.connections.values():
                    connection.flow_rate = connection.calc_flow_rate()
                return self

            if n_iter >= max_iter:
                raise RuntimeError(f'Simulation not converged in {max_iter} iteration(s)!')
