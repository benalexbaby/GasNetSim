#   #!/usr/bin/env python
#   -*- coding: utf-8 -*-
#   ******************************************************************************
#     Copyright (c) 2023.
#     Developed by Yifei Lu
#     Last change on 1/5/23, 9:16 AM
#     Last change by yifei
#    *****************************************************************************
import pandas as pd
import numpy as np
from tqdm import tqdm
from numpy.testing import assert_almost_equal

from GasNetSim.components import *
from GasNetSim.components.utils.pipeline_function.outlet_temperature import *

from scipy.constants import bar, R, atm

import pytest


def test_pipeline_outlet_temperature_calculation():

    # test_nodes = {1: Node(node_index=1, pressure_pa=50 * bar, temperature=300, node_type='reference'),
    #               2: Node(node_index=2, pressure_pa=48 * bar, temperature=300, node_type='reference')}
    # test_pipeline = Pipeline(test_nodes[1], test_nodes[2], diameter=0.5, length=500*1e3)

    gas_mixture = Mixture(T=300, P=50 * bar, zs={"methane": 0.9, "hydrogen": 0.1})
    heat_transfer_coefficient = 3.69 * 10
    mass_flow_rate = 55  # kg/s
    Cp = gas_mixture.Cp
    diameter = 0.5
    JT = gas_mixture.JT
    Z = gas_mixture.Z
    R_specific = R / gas_mixture.MW
    friction_factor = 0.01
    average_pressure = 50*bar
    beta = calculate_beta_coefficient(heat_transfer_coefficient, mass_flow_rate, Cp, diameter)
    gamma = calculate_gamma_coefficient(JT, Z, R_specific, friction_factor, mass_flow_rate, average_pressure, diameter)
    T_outlet = calculate_pipeline_outlet_temperature(beta, gamma, 288.15, 150000, 300)

    dx = 1000
    position = np.arange(0, 150000 + dx, dx)
    Xs_EF = np.zeros(len(position))  # Euler forward
    Xs_EB = np.zeros(len(position))  # Euler backward
    Xs_TR = np.zeros(len(position))  # Trapizoid rule

    Xs_EF[0] = Xs_EB[0] = Xs_TR[0] = 300

    a = gamma - beta
    b = beta * 288.15
    for i in range(1, 151):
        Xs_EF[i] = Xs_EF[i - 1] + dx * (a * Xs_EF[i - 1] + b)
        Xs_EB[i] = 1 / (1 - dx * a) * (Xs_EB[i - 1] + dx * b)
        Xs_TR[i] = (1 + dx * a / 2) / (1 - dx * a / 2) * Xs_TR[i - 1] + (dx * b) / (1 - dx * a / 2)

    assert_almost_equal(T_outlet, Xs_EF[-1])
    assert_almost_equal(T_outlet, Xs_EB[-1])
    assert_almost_equal(T_outlet, Xs_TR[-1])


if __name__ == "__main__":
    pipeline_temperature_profile = list()

    network_nodes = {1: Node(node_index=1, pressure_pa=50 * bar, temperature=300, node_type='reference'),
                     2: Node(node_index=2, flow=20)}

    for length in tqdm(range(100)):

        network_pipes = {1: Pipeline(network_nodes[1], network_nodes[2], diameter=0.5, length=(length + 1) * 1e3)}
        gas_network = Network(network_nodes, network_pipes)
        try:
            gas_network.simulation()
        except RuntimeError:
            logging.error(f"Simulation does not converge for pipeline length {length + 1} km!")
        t = network_pipes[1].calc_pipe_outlet_temp()
        pipeline_temperature_profile.append(t)

    plt.figure()
    plt.plot(pd.Series(pipeline_temperature_profile))
    plt.xlabel("Pipeline length [km]")
    plt.ylabel("Temperature [K]")
    plt.title("Temperature profile alongside the pipeline")
    plt.show()
