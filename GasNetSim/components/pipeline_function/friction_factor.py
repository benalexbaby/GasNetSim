#  #!/usr/bin/env python
#  -*- coding: utf-8 -*-
#  ******************************************************************************
#    Copyright (c) 2022.
#    Developed by Yifei Lu
#    Last change on 1/17/22, 11:17 AM
#    Last change by yifei
#   *****************************************************************************
import numpy as np
import math
from scipy.optimize import fsolve
import warnings


def reynold_number(diameter, velocity, rho, viscosity):
    """

    :param diameter: pipe diameter (m)
    :param velocity: fluid velocity (m/s)
    :param rho: fluid density (kg/m3)
    :param viscosity: fluid viscosity(kg/(m*s))
    :return: Reynold number
    """
    return diameter * velocity * rho / viscosity


def reynold_number_simp(diameter, sg, q, viscosity):
    """
    A simplified method to calculate the Reynolds number based on the volumetric flow rate
    :param diameter: pipe diameter (m)
    :param sg: gas specific gravity
    :param q: gas flow rate (sm3/s)
    :param viscosity: fluid viscosity (Pa*s)
    :return:
    """
    pb = 101.325  # kPa
    tb = 288.15  # K
    return 49.44 * q * sg * pb / (viscosity * diameter * tb) * (24*3600)


def hagen_poiseuille(N_re):
    """
    Friction factor in Laminar zone can be calculated using Hagen-Poiseuille method
    :param N_re:
    :return:
    """
    if N_re >= 2100:
        warnings.warn("You are using Hagen-Poiseuille friction model for a non-laminar flow!")
    return 64 / N_re


def nikuradse(d, epsilon):
    d *= 1000
    return 1 / (2 * math.log(d / epsilon, 10) + 1.14) ** 2


def von_karman_prandtl(N_re):
    """

    :param N_re: Reynolds number
    :return: von Karman - Prandtl friction factor
    """
    def func(f): return 2 * math.log(N_re * math.sqrt(f), 10) - 0.8 - 1 / math.sqrt(f)
    f_init_guess = np.array(0.01)
    friction_factor = fsolve(func, f_init_guess)
    return friction_factor


def colebrook_white(epsilon, d, N_re):
    """

    :param epsilon:
    :param d:
    :param N_re:
    :return:
    """
    d *= 1000
    def func(f): return -2 * math.log(epsilon/d/3.71 + 2.51/N_re/math.sqrt(f), 10) - 1 / math.sqrt(f)
    f_init_guess = np.array(0.01)
    friction_factor = fsolve(func, f_init_guess)
    return friction_factor


def colebrook_white_hofer_approximation(N_re, d, epsilon):
    return (-2 * math.log(4.518/N_re * math.log(N_re/7, 10) + epsilon/3.71/d, 10))**(-2)


def nikuradse_from_CWH(epsilon, d):
    return (-2 * math.log(epsilon/3.71/d)) ** (-2)


def chen(epsilon, d, N_re):
    d *= 1000
    _ = epsilon/d/3.7065 - 5.0452/N_re * math.log((((epsilon/d)**1.1096)/2.8257 + (7.149/N_re)**0.8961), 10)
    return 1/(4 * math.log(_, 10) ** 2)


def weymouth(d):
    return 0.0093902 / (d ** (1 / 3))


if __name__ == "__main__":
    from thermo import Mixture
    from collections import OrderedDict
    import math
    import matplotlib.pyplot as plt

    gas_comp = OrderedDict([('methane', 0.96522),
                            ('nitrogen', 0.00259),
                            ('carbon dioxide', 0.00596),
                            ('ethane', 0.01819),
                            ('propane', 0.0046),
                            ('isobutane', 0.00098),
                            ('butane', 0.00101),
                            ('2-methylbutane', 0.00047),
                            ('pentane', 0.00032),
                            ('hexane', 0.00066)])

    Nre_res = list()
    Nre_res_simp = list()

    for p in range(1, 100):
        gas_mixture = Mixture(T=288.15, P=p*101325, zs=gas_comp)

        gas_mix_viscosity = gas_mixture.mu
        gas_mix_density = gas_mixture.rho
        pipe_diameter = 0.76  # m
        volumetric_flow_rate = 20  # m3/s
        flow_velocity = volumetric_flow_rate / (math.pi * (pipe_diameter/2)**2) * p # m/s
        gas_mix_specific_gravity = gas_mixture.SG

        Nre = reynold_number(pipe_diameter, flow_velocity, gas_mix_density, gas_mix_viscosity)
        Nre_simple = reynold_number_simp(pipe_diameter*1000, gas_mix_specific_gravity, volumetric_flow_rate, gas_mix_viscosity)
        # print("Reynolds' number is: {}".format(Nre))
        # print("Reynolds' number calculated with simplified method is: {}".format(Nre_simple))
        Nre_res.append(Nre)
        Nre_res_simp.append(Nre_simple)

    plt.figure()
    plt.plot(Nre_res)
    plt.plot(Nre_res_simp)
    plt.show()