#  #!/usr/bin/env python
#  -*- coding: utf-8 -*-
#  ******************************************************************************
#    Copyright (c) 2022.
#    Developed by Yifei Lu
#    Last change on 3/9/22, 11:18 AM
#    Last change by yifei
#   *****************************************************************************
import math


def calculate_beta_coefficient(ul, qm, cp, D):
    """

    :param ul: Heat transfer coefficient [W m^(-2) K^(-1)]
    :param qm: Mass flow [kg/s]
    :param cp: Specific heat capacity at constant pressure [J/(kg K)]
    :param D: pipeline diameter [m]
    :return:
    """
    return ul/(qm*cp) * math.pi * D


def calculate_gamma_coefficient(mu_jt, Z, R_specific, f, qm, p_avg, D):
    """

    :param mu_jt: Joule–Thomson coefficient [K/Pa]
    :param Z: Compressibility factor [-]
    :param R_specific:
    :param f: friction factor [-]
    :param qm: Mass flow rate [kg/s]
    :param p_avg: Average pressure [Pa]
    :param D: Pipeline diameter [m]
    :return: Gamma coefficient used for the calculation of pipeline outlet temperature
    """
    A = math.pi * (D/2)**2
    return mu_jt * Z * R_specific * f * qm * abs(qm) / (2*p_avg*D*A**2)


def calculate_pipeline_outlet_temperature(beta, gamma, Ts, L, T1):
    """

    :param beta: Beta coefficient used for the calculation of pipeline outlet temperature
    :param gamma: Gamma coefficient used for the calculation of pipeline outlet temperature
    :param Ts: Surrounding temperature [K]
    :param L: Pipeline length [m]
    :param T1: Temperature at pipeline inlet [K]
    :return: Temperature at pipeline outlet [K]
    """
    return beta/(beta-gamma)*(Ts - Ts * math.exp((gamma-beta)*L)) + T1*math.exp((gamma-beta)*L)


def calculate_steady_state_outlet_temperature(beta, gamma, Ts):
    return beta * Ts / (beta - gamma)


if __name__ == "__main__":
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.constants import bar, R

    # from GasNetSim.components.utils.gas_mixture.thermo.thermo import Mixture
    from GasNetSim.components.utils.gas_mixture.gas_mixture import GasMixture

    gas_mixture = GasMixture(temperature=300,
                             pressure=50 * bar,
                             composition={"methane": 0.9, "hydrogen": 0.1})
    heat_transfer_coefficient = 3.69 * 10
    mass_flow_rate = 55  # kg/s
    Cp = gas_mixture.heat_capacity_constant_pressure
    diameter = 0.5
    JT = gas_mixture.joule_thomson_coefficient
    Z = gas_mixture.compressibility
    R_specific = gas_mixture.R_specific
    friction_factor = 0.01
    average_pressure = 50*bar
    beta = calculate_beta_coefficient(heat_transfer_coefficient, mass_flow_rate, Cp, diameter)
    gamma = calculate_gamma_coefficient(JT, Z, R_specific, friction_factor, mass_flow_rate, average_pressure, diameter)
    T_outlet = calculate_pipeline_outlet_temperature(beta, gamma, 288.15, 150000, 300)
    print(beta, gamma)
    print(T_outlet)

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

    plt.figure()
    plt.plot(Xs_EF, label="Euler forward")
    plt.plot(Xs_EB, label="Euler backward")
    plt.plot(Xs_TR, label="Trapizoid rule")
    plt.xlabel("Distance [km]")
    plt.ylabel("Temperature [K]")
    plt.legend()
    plt.show()
