#  #!/usr/bin/env python
#  -*- coding: utf-8 -*-
#  ******************************************************************************
#    Copyright (c) 2022.
#    Developed by Yifei Lu
#    Last change on 3/14/22, 9:41 PM
#    Last change by yifei
#   *****************************************************************************
import pandas as pd
import numpy as np
from collections import OrderedDict
from pathlib import Path

from ..components.pipeline import Pipeline
from ..components.node import Node
from ..components.network import Network
from ..utils.exception import *


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


def read_resistances(path_to_file: Path) -> dict:
    """

    :param path_to_file:
    :return:
    """
    resistances = dict()
    return resistances


def create_network_from_csv(path_to_folder: Path) -> Network:
    """

    :param path_to_folder:
    :return:
    """
    all_files = list(path_to_folder.glob('*.csv'))
    nodes = read_nodes(Path('./' + '_'.join(all_files[0].stem.split('_')[:-1]) + '_nodes.csv'))
    for file in all_files:
        file_name = file.stem
        if 'node' in file_name:
            pass
        elif 'pipeline' in file_name:
            pipelines = read_pipelines(file, nodes)
        elif 'compressor' in file_name:
            compressors = read_compressors(file)
        elif 'resistance' in file_name:
            resistances = read_resistances(file)
        else:
            raise FileNameError(f'Please check the file name {file_name}.csv')

    return Network(nodes=nodes, pipelines=pipelines)
