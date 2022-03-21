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

from .node import *
from .pipeline import *
from .utils import *

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class Network:
    """
    Network class
    """

    def __init__(self, nodes: dict, pipelines: dict, compressors=None, resistances=None):
        """

        :param nodes:
        :param pipelines:
        """
        self.nodes = nodes
        self.pipelines = pipelines
        self.compressors = compressors
        self.resistances = resistances
        self.p_ref_nodes = [x for x in self.nodes.values() if x.node_type == 'reference']
        self.t_ref_nodes = [x for x in self.nodes.values() if x.node_type == 'reference']
        self.p_ref_nodes_index = [x for x, y in self.nodes.items() if y.node_type == 'reference']
        self.t_ref_nodes_index = [x for x, y in self.nodes.items() if y.node_type == 'reference']

    def find_reference_nodes(self) -> Tuple[list, list]:
        """
        Find reference nodes, where pressure are pre-set
        :return: List of Node class of pressure-referenced nodes and List of Node class of temperature-referenced nodes
        """
        pressure_ref_nodes = list()
        temperature_ref_nodes = list()

        for node in self.nodes.values():
            if node.pressure is not None:
                pressure_ref_nodes.append(node)
            if node.temperature is not None:
                temperature_ref_nodes.append(node)

        return pressure_ref_nodes, temperature_ref_nodes

    def find_ref_nodes_index(self):
        """
        To get indices of the referenced nodes
        :return: List of pressure-referenced nodes indices and list of temperature-referenced nodes indices
        """
        pressure_ref_nodes_index = list()
        temperature_ref_nodes_index = list()

        for node_index, node in self.nodes.items():
            if node.pressure is not None:
                pressure_ref_nodes_index.append(node_index)
            if node.temperature is not None:
                temperature_ref_nodes_index.append(node_index)

        return pressure_ref_nodes_index, temperature_ref_nodes_index

    @property
    def demand_nodes_supply_pipelines(self):
        """

        :return:
        """
        nodal_supply_pipelines = dict()
        for i_pipe, pipe in self.pipelines.items():
            if nodal_supply_pipelines.get(pipe.outlet_index) is not None:
                nodal_supply_pipelines[pipe.outlet_index].append(i_pipe)
            else:
                nodal_supply_pipelines[pipe.outlet_index] = [i_pipe]
        return OrderedDict(sorted(nodal_supply_pipelines.items()))

    def convert_energy_flow_to_volumetric_flow(self, base='HHV'):
        for node in self.nodes.values():
            gas_comp = node.get_mole_fraction()
            standard_density = Mixture(P=101325, T=288.15, zs=gas_comp).rho
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

    def connection_matrix(self, sparse_matrix=False):
        # TODO change the index number
        n_nodes = len(self.nodes.values())
        pipelines = self.pipelines
        compressors = self.compressors
        resistances = self.resistances
        if sparse_matrix:
            row_ind = list()
            col_ind = list()
            data = list()

        # Build a matrix to show the connection between nodes
        connection = np.zeros((n_nodes, n_nodes))
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

        return connection

    def pressure_initialization(self):
        nodes = self.nodes
        pipelines = self.pipelines
        max_length = max([x.length for x in pipelines.values()])
        max_flow = max([x.flow for x in nodes.values() if x.flow is not None])
        pressure_init = [node.pressure for node in nodes.values()]
        # pipeline_with_missing_pressure = copy.deepcopy(pipelines)
        pressure_init_old = list()

        while pressure_init != pressure_init_old:
            pressure_init_old = copy.deepcopy(pressure_init)
            # pipeline_initialized = list()
            for i_pipe, pipe in pipelines.items():
                i = pipe.inlet_index - 1
                j = pipe.outlet_index - 1
                length = pipe.length
                flow = pipe.outlet.flow
                if pressure_init[i] is None and pressure_init[j] is None:
                    pass
                elif pressure_init[j] is None or pressure_init[i] == pressure_init[j]:
                    if length/max_length < 0.01:
                        pressure_init[j] = pressure_init[i] * 0.999999
                    else:
                        # pressure_init[j] = pressure_init[i] * (1 - 0.05 * (length/max_length)**0.5 * (flow/max_flow))
                        pressure_init[j] = pressure_init[i] * 0.98
                elif pressure_init[j] is not None and pressure_init[i] is not None:
                    if length / max_length < 0.01:
                        pressure_init[j] = min(pressure_init[j], pressure_init[i] * 0.99999)
                    else:
                        # pressure_init[j] = min(pressure_init[j],
                        #                        pressure_init[i] * (1 - 0.05 * (length/max_length)**0.5 * (flow/max_flow)))
                        pressure_init[j] = min(pressure_init[j], pressure_init[i] * 0.98)
                elif pressure_init[i] is None and pressure_init[j] is not None:
                    if length/max_length < 0.01:
                        pressure_init[i] = pressure_init[j] / 0.99999
                    else:
                        # pressure_init[i] = pressure_init[j] / (1 - 0.05 * (length/max_length)**0.5 * (flow /max_flow))
                        pressure_init[i] = pressure_init[j] / 0.98

        return pressure_init

    def newton_raphson_initialization(self):
        """
        Initialization for NR-solver, where the nodal pressures are initialized as 0.98 of inlet pressure and nodal
        temperatures are the same as pipe surrounding temperatures
        :return: Network initial conditions for NR-solver
        """

        nodes = self.nodes
        pipelines = self.pipelines

        n_nodes = len(nodes)

        # Build a matrix to show the connection between nodes
        connection = np.zeros((n_nodes, n_nodes))
        for pipe in pipelines.values():
            i = pipe.inlet_index - 1
            j = pipe.outlet_index - 1
            connection[i][j] = 1
            connection[j][i] = 1

        # TODO consider the case where ref_nodes do not start with index 0
        p_ref_nodes = self.p_ref_nodes_index

        for node in self.nodes.values():
            HHV = calc_heating_value(node.gas_mixture)
            if node.flow_type == 'volumetric':
                pass
            elif node.flow_type == 'energy':
                gas_comp = node.get_mole_fraction()
                node.flow = node.flow / HHV * 1e6 / Mixture(zs=gas_comp, T=288.15, P=101325).rho
                print(node.flow)
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
        # logger.debug(pressure_init)
        # for i in range(len(pressure_init)):
        #     for j in range(len(pressure_init)):
        #         if connection[i][j] != 0 and pressure_init[i] is None and pressure_init[j] is None:
        #             pass
        #         elif connection[i][j] != 0 and (pressure_init[j] is None or pressure_init[i] == pressure_init[j]):
        #             pressure_init[j] = pressure_init[i] * 0.98
        #         elif connection[i][j] != 0 and pressure_init[i] is None:
        #             pressure_init[i] = pressure_init[j] / 0.98

        for i in range(len(nodal_flow_init)):
            # TODO change to number of non-reference nodes
            nodes[i + 1].pressure = pressure_init[i]
            nodes[i + 1].flow = nodal_flow_init[i]
            nodes[i + 1].volumetric_flow = nodal_flow_init[i]
            nodes[i + 1].convert_volumetric_to_energy_flow()
            nodes[i + 1].temperature = temperature_init[i]

        for index, pipe in pipelines.items():
            pipe.inlet = nodes[pipe.inlet_index]
            pipe.outlet = nodes[pipe.outlet_index]

            # try:
            #     pipe.update_gas_mixture()
            # except TypeError:
            #     raise ValueError("Something is wrong with gas mixture initialization...")

        return nodal_flow_init, pressure_init, temperature_init

    def jacobian_matrix(self):

        pipelines = self.pipelines
        nodes = self.nodes
        p_ref_nodes, t_ref_nodes = self.p_ref_nodes, self.t_ref_nodes
        # p_ref_nodes_index, t_ref_nodes_index = self.p_ref_nodes_index, self.t_ref_nodes_index
        n_nodes = len(nodes)
        n_p_nodes = n_nodes - len(p_ref_nodes)

        n_t_nodes = n_nodes - len(t_ref_nodes)

        n_non_ref_nodes = n_nodes - len(p_ref_nodes)

        jacobian_mat = np.zeros((n_non_ref_nodes, n_non_ref_nodes), dtype=np.float)
        flow_mat = np.zeros((n_nodes, n_nodes), dtype=np.float)

        for pipe in pipelines.values():
            flow_mat[pipe.inlet_index - 1][pipe.outlet_index - 1] = - pipe.calc_flow_rate()
            flow_mat[pipe.outlet_index - 1][pipe.inlet_index - 1] = pipe.calc_flow_rate()

            if pipe.inlet_index in p_ref_nodes and pipe.outlet_index in p_ref_nodes:
                pass
            else:
                i = pipe.inlet_index - 1 - len(p_ref_nodes)
                j = pipe.outlet_index - 1 - len(p_ref_nodes)

                slope_corr = pipe.calc_pipe_slope_correction()
                p1 = pipe.inlet.pressure
                p2 = pipe.outlet.pressure
                phy_char = pipe.calc_physical_char_gas_pipe()
                temp_var = (abs(p1 ** 2 - p2 ** 2 - slope_corr)) ** (-0.5)

                if i >= 0 and j >= 0:
                    jacobian_mat[i][j] = phy_char * p2 * temp_var * 2
                    jacobian_mat[j][i] = phy_char * p1 * temp_var * 2

                if i >= 0:
                    jacobian_mat[i][i] += - phy_char * p1 * temp_var * 2
                if j >= 0:
                    jacobian_mat[j][j] += - phy_char * p2 * temp_var * 2

        return jacobian_mat, flow_mat

    def simulation(self):
        print([x.flow for x in self.nodes.values()])
        ref_nodes = self.p_ref_nodes_index

        n_nodes = len(self.nodes.keys())
        connection_matrix = np.zeros((n_nodes, n_nodes))
        for pipe in self.pipelines.values():
            i = pipe.inlet_index - 1
            j = pipe.outlet_index - 1
            connection_matrix[i][j] = 1
            connection_matrix[j][i] = 1

        init_f, init_p, init_t = self.newton_raphson_initialization()

        max_iter = 100
        n_iter = 0
        n_non_ref_nodes = n_nodes - len(ref_nodes)

        f = np.array(init_f)
        p = np.array(init_p)
        t = np.array(init_t)
        print(f'Initial pressure: {p}')
        # print(f'Initial flow: {f}')

        for i in range(len(init_f)):
            # TODO change to number of non-reference nodes
            self.nodes[i + 1].pressure = init_p[i]
            self.nodes[i + 1].flow = init_f[i]
            self.nodes[i + 1].temperature = init_t[i]
            self.nodes[i + 1].update_gas_mixture()

        for index, pipe in self.pipelines.items():
            pipe.inlet = self.nodes[pipe.inlet_index]
            pipe.outlet = self.nodes[pipe.outlet_index]
            pipe.update_gas_mixture()

        delta_flow = 0

        record = list()

        while n_iter <= max_iter:

            # Calculate nodal inflow gas mixture composition
            nodal_gas_inflow_comp = dict()
            demand_node_supply_pipelines = self.demand_nodes_supply_pipelines
            for i_node, node in self.nodes.items():
                if i_node in ref_nodes:
                    nodal_gas_inflow_comp[i_node] = node.get_mole_fraction()
                elif self.nodes[i_node].flow < 0:  # supply nodes
                    nodal_gas_inflow_comp[i_node] = node.get_mole_fraction()
                else:
                    nodal_gas_inflow_comp[i_node] = 0

                    supply_pipelines = [self.pipelines[j] for j in demand_node_supply_pipelines[i_node]]

                    total_inflow_comp = dict()
                    total_inflow = 0
                    total_inflow_temperature_times_flow_rate = 0
                    for pipe in supply_pipelines:
                        pipe.gas_mixture = pipe.inlet.gas_mixture
                        inflow_rate = pipe.calc_flow_rate()
                        if inflow_rate > 0:
                            # Sum up flow rate * temperature
                            total_inflow_temperature_times_flow_rate += inflow_rate * pipe.calc_pipe_outlet_temp()
                            total_inflow += inflow_rate

                            gas_comp = pipe.get_mole_fraction()
                            # create a OrderedDict to store gas flow fractions
                            gas_flow_comp = OrderedDict({gas: comp * inflow_rate for gas, comp in gas_comp.items()})
                            for gas, comp in gas_flow_comp.items():
                                if total_inflow_comp.get(gas) is None:
                                    total_inflow_comp[gas] = comp
                                else:
                                    total_inflow_comp[gas] += comp
                            # if nodal_gas_inflow_comp[i_node] == 0:
                            #     nodal_gas_inflow_comp[i_node] = gas_flow_comp
                            # else:
                            #     nodal_gas_inflow_comp[i_node] = dict(Counter(nodal_gas_inflow_comp[i_node]) + Counter(gas_flow_comp))
                    total = sum(total_inflow_comp.values(), 0.0)
                    nodal_gas_inflow_comp[i_node] = {k: v / total for k, v in total_inflow_comp.items()}
                    if total_inflow != 0:
                        temperature_new = total_inflow_temperature_times_flow_rate / total_inflow
                        print(temperature_new)

                        # update node temperature
                        self.nodes[i_node].temperature = temperature_new

                    try:
                        self.nodes[i_node].gas_mixture = Mixture(zs=nodal_gas_inflow_comp[i_node],
                                                                 T=self.nodes[i_node].temperature,
                                                                 P=self.nodes[i_node].pressure)
                    except Exception:
                        print(i_node)
                        print(i_pipe)
                        print(nodal_gas_inflow_comp[i_node])

            j_mat, f_mat = self.jacobian_matrix()
            j_mat_inv = np.linalg.inv(j_mat)

            delta_flow = f - np.dot(f_mat, np.ones(n_nodes))

            delta_flow = [delta_flow[i] for i in range(len(delta_flow)) if i + 1 not in ref_nodes]

            # Update volumetric flow rate target
            for n in self.nodes.values():
                n.convert_energy_to_volumetric_flow()
            f = [x.volumetric_flow if x.flow is not None else 0 for x in self.nodes.values()]

            try:
                delta_p = np.dot(j_mat_inv, delta_flow)
            except ValueError:
                print(delta_flow)

            p += np.concatenate((np.array([0] * len(ref_nodes)), delta_p), axis=None)

            for i in self.nodes.keys():
                if i not in ref_nodes:
                    self.nodes[i].pressure = p[i - 1]

            for i_pipe, pipe in self.pipelines.items():
                pipe.inlet = self.nodes[pipe.inlet_index]
                pipe.outlet = self.nodes[pipe.outlet_index]

            record.append(delta_p)

            n_iter += 1
            j_mat, f_mat = self.jacobian_matrix()
            delta_flow = f[len(ref_nodes):] - np.dot(f_mat, np.ones(n_nodes))[len(ref_nodes):]

            print(max([abs(x) for x in (delta_flow/f[len(ref_nodes):])]))
            print(delta_p)

            if max([abs(x) for x in (delta_flow/f[len(ref_nodes):])]) <= 0.001:
                logger.info(f'Simulation converges in {n_iter} iterations.')
                logger.info(p)
                logger.debug(init_p)
                pipe_h2_fraction = list()
                for i_pipe, pipe in self.pipelines.items():
                    logger.debug(f'Pipeline index: {i_pipe}')
                    logger.debug(f'Pipeline flow rate: {pipe.flow_rate}')
                    logger.debug(f'Gas mixture composition: {pipe.get_mole_fraction()}')
                    try:
                        pipe_h2_fraction.append(pipe.get_mole_fraction()['hydrogen'] * 100)
                    except KeyError:
                        pipe_h2_fraction.append(0)
                print(pipe_h2_fraction)
                for pipe in self.pipelines.values():
                    pipe.flow_rate = pipe.calc_flow_rate()
                return self

            if n_iter >= max_iter:
                raise RuntimeError(f'Simulation not converged in {max_iter} iteration(s)!')
