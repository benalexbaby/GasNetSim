#   #!/usr/bin/env python
#   -*- coding: utf-8 -*-
#   ******************************************************************************
#     Copyright (c) 2023.
#     Developed by Yifei Lu
#     Last change on 5/31/23, 4:14 PM
#     Last change by yifei
#    *****************************************************************************
import pandas as pd
import numpy as np
from collections import OrderedDict
from pathlib import Path

from GasNetSim.components.pipeline import Pipeline, Resistance, ShortPipe
from GasNetSim.components.node import Node
from GasNetSim.components.network import Network
from GasNetSim.utils.exception import *


def read_nodes(path_to_file: Path) -> dict:
    """

    :param path_to_file:
    :return:
    """
    nodes = dict()
    df_node = pd.read_csv(path_to_file, delimiter=';')
    df_node = df_node.replace({np.nan: None})

    for row_index, row in df_node.iterrows():
        if row['gas_composition'] is not None:
            row['gas_composition'] = OrderedDict(eval(row['gas_composition']))
        nodes[row['node_index']] = Node(node_index=row['node_index'],
                                        pressure_pa=row['pressure_pa'],
                                        flow=row['flow_sm3_per_s'],
                                        temperature=row['temperature_k'],
                                        altitude=row['altitude_m'],
                                        gas_composition=row['gas_composition'],
                                        node_type=row['node_type'],
                                        flow_type=row['flow_type'])
    return nodes


def read_pipelines(path_to_file: Path, network_nodes: dict) -> dict:
    """

    :param path_to_file:
    :param network_nodes:
    :return:
    """
    pipelines = dict()
    df_pipe = pd.read_csv(path_to_file, delimiter=';')
    df_pipe = df_pipe.replace({np.nan: None})

    for row_index, row in df_pipe.iterrows():
        pipelines[row['pipeline_index']] = Pipeline(inlet=network_nodes[row['inlet_index']],
                                                    outlet=network_nodes[row['outlet_index']],
                                                    diameter=row['diameter_m'],
                                                    length=row['length_m'])
    return pipelines


def read_compressors(path_to_file: Path) -> dict:
    """

    :param path_to_file:
    :return:
    """
    compressors = dict()
    return compressors


def read_resistances(path_to_file: Path, network_nodes: dict) -> dict:
    """

    :param path_to_file:
    :param network_nodes:
    :return:
    """
    resistances = dict()
    df_resistance = pd.read_csv(path_to_file, delimiter=';')
    df_resistance = df_resistance.replace({np.nan: None})

    for row_index, row in df_resistance.iterrows():
        resistances[row['resistance_index']] = Resistance(inlet=network_nodes[row['inlet_index']],
                                                          outlet=network_nodes[row['outlet_index']],
                                                          resistance=row['resistance'])
    return resistances


def read_shortpipes(path_to_file: Path, network_nodes: dict) -> dict:
    """

    :param path_to_file:
    :param network_nodes:
    :return:
    """
    shortpipes = dict()
    df_shortpipes = pd.read_csv(path_to_file, delimiter=';')
    df_shortpipes = df_shortpipes.replace({np.nan: None})

    for row_index, row in df_shortpipes.iterrows():
        shortpipes[row['shortpipe_index']] = ShortPipe(inlet=network_nodes[row['inlet_index']],
                                                       outlet=network_nodes[row['outlet_index']])
    return shortpipes


def create_network_from_csv(path_to_folder: Path) -> Network:
    """

    :param path_to_folder:
    :return:
    """
    all_files = list(path_to_folder.glob('*.csv'))
    nodes = read_nodes(Path('./' + '_'.join(all_files[0].stem.split('_')[:-1]) + '_nodes.csv'))

    network_components = {'node': nodes,  # the dataset should have at least node
                          'pipeline': None,
                          'compressor': None,
                          'resistance': None,
                          'shortpipe': None}

    for file in all_files:
        file_name = file.stem
        if 'node' in file_name:
            pass
        elif 'pipeline' in file_name:
            pipelines = read_pipelines(file, nodes)
            network_components['pipeline'] = pipelines
        elif 'compressor' in file_name:
            compressors = read_compressors(file)
            network_components['compressor'] = compressors
        elif 'resistance' in file_name:
            resistances = read_resistances(file, nodes)
            network_components['resistance'] = resistances
        elif 'shortpipe' in file_name:
            shortpipes = read_shortpipes(file, nodes)
            network_components['shortpipe'] = shortpipes
        else:
            raise FileNameError(f'Please check the file name {file_name}.csv')

    return Network(nodes=network_components['node'],
                   pipelines=network_components['pipeline'],
                   compressors=network_components['compressor'],
                   resistances=network_components['resistance'],
                   shortpipes=network_components['shortpipe'])
