#   #!/usr/bin/env python
#   -*- coding: utf-8 -*-
#   ******************************************************************************
#     Copyright (c) 2022.
#     Developed by Yifei Lu
#     Last change on 7/28/22, 4:19 PM
#     Last change by yifei
#    *****************************************************************************
import pandas as pd
from pathlib import Path
import logging
import copy
from tqdm import tqdm

from ..components.network import Network
from ..components.utils.utils import plot_network_demand_distribution


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
logger.setLevel(level=logging.WARNING)


def read_profiles(file, sep=";"):
    profiles = pd.read_csv(Path(file), sep=sep)
    logger.info(f'Reading profiles from {file}.')
    return profiles


def print_progress_bar(iteration, total, prefix='', suffix='', decimals=1, length=100, fill='█'):
    """
    Call in a loop to create terminal progress bar.
    the code is mentioned in : https://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filled_length = int(length * iteration // total)
    bar = fill * filled_length + '-' * (length - filled_length)
    # logger.info('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix))
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end="")
    # Print New Line on Complete
    if iteration == total:
        print("\n")


def check_profiles(profiles):
    if profiles['time'].dtype == int:
        print()
    elif profiles['time'].dtype == pd.Timestamp:
        print()
        # df['time'] = df['time'].apply(lambda x: pd.Timestamp.now() + pd.Timedelta(seconds=x))


def run_snapshot(network, tol=0.01):
    # plot_network_demand_distribution(network)
    network = network.simulation(tol=tol)
    return network


def not_converged(time_step, ts_variables):
    logger.error(f'CalculationNotConverged at time step {time_step}.')
    if not ts_variables["continue_on_divergence"]:
        raise ts_variables['errors'][0]


def update_network_topology(network):
    full_network = copy.deepcopy(network)  # make a copy of the fully connected network
    network_nodes = full_network.nodes
    network_pipes = full_network.pipelines
    network_resistances = full_network.resistances

    # Check which nodes need to be removed
    removed_nodes = dict()
    remaining_nodes = dict()
    for i, node in list(network_nodes.items()):
        if node.flow == 0:
            removed_nodes[i] = node
        else:
            remaining_nodes[i - len(removed_nodes)] = node

    # Check which pipelines need to be removed
    removed_pipes = dict()
    remaining_pipes = dict()
    for i, pipe in list(network_resistances.items()):
        if (pipe.inlet_index in removed_nodes.keys()) or (pipe.outlet_index in removed_nodes.keys()):
            pipe.valve = 1
        else:
            pipe.valve = 0
        if pipe.valve == 1:
            removed_pipes[i] = pipe
        else:
            pipe.inlet_index = list(remaining_nodes.keys())[list(remaining_nodes.values()).index(pipe.inlet)]
            pipe.outlet_index = list(remaining_nodes.keys())[list(remaining_nodes.values()).index(pipe.outlet)]
            remaining_pipes[i - len(removed_pipes)] = pipe

    return Network(nodes=remaining_nodes, pipelines=None, resistances=remaining_pipes)


def run_time_series(network, file=None, sep=";", profile_type="energy"):
    # create a copy of the input network
    full_network = copy.deepcopy(network)
    results_to_save = ["nodal_pressure", "pipeline_flowrate", "nodal_gas_composition"]
    results = dict([(k, []) for k in results_to_save])

    # read profile
    if file is not None:
        profiles = read_profiles(file, sep=sep)
        time_steps = profiles.index
    else:
        time_steps = range(5)  # test with 5 fictitious time steps

    # create error log to record the time step indices where error occurs
    error_log = list()

    pressure_prev = None

    for t in tqdm(time_steps):
        full_network = copy.deepcopy(network)
        full_network.pressure_prev = pressure_prev  # Nodal pressure values at previous time step
        if pressure_prev is None:
            full_network.run_initialization = True
        else:
            full_network.run_initialization = False
            full_network.assign_pressure_values(pressure_prev)

        for i in full_network.nodes.keys():
            if i in full_network.reference_nodes:
                full_network.nodes[i].volumetric_flow = None
                full_network.nodes[i].energy_flow = None
            else:
                try:
                    if profile_type == "volumetric":
                        full_network.nodes[i].volumetric_flow = profiles[str(i)][t]
                        full_network.nodes[i].convert_volumetric_to_energy_flow()
                    elif profile_type == "energy":
                        full_network.nodes[i].energy_flow = profiles[str(i)][t]
                        full_network.nodes[i].convert_energy_to_volumetric_flow()
                    else:
                        raise ValueError(f"Unknown profile type {profile_type}!")
                    # full_network.nodes[i].demand_type = 'energy'
                except KeyError:
                    print(f"Node index {i} is not found!")
        # simplified_network = update_network_topology(full_network)
        try:
            # network = run_snapshot(simplified_network)
            for n in full_network.nodes.values():
                if n.volumetric_flow is not None and n.volumetric_flow < 0:
                    print(n.volumetric_flow)
            full_network = copy.deepcopy(run_snapshot(full_network))
            pressure_prev = full_network.save_pressure_values()
        except RuntimeError:
            # error_log.append([simplified_network, profiles.iloc[t]])
            error_log.append([full_network, profiles.iloc[t]])

        results = save_time_series_results(full_network, results, nodal_pressure=True, pipeline_flowrate=True,
                                           nodal_gas_composition=True)

    return results


def save_time_series_results(network, results, nodal_pressure=True, pipeline_flowrate=True, nodal_gas_composition=True):
    if nodal_pressure:
        results["nodal_pressure"].append([node.pressure for node in network.nodes.values()])

    if nodal_gas_composition:
        results["nodal_gas_composition"].append([node.gas_mixture.composition for node in network.nodes.values()])

    if pipeline_flowrate:
        results["pipeline_flowrate"].append([pipe.flow_rate for pipe in network.pipelines.values()])

    return results