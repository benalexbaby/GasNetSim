#  #!/usr/bin/env python
#  -*- coding: utf-8 -*-
#  ******************************************************************************
#    Copyright (c) 2022.
#    Developed by Yifei Lu
#    Last change on 3/21/22, 3:22 PM
#    Last change by yifei
#   *****************************************************************************
import math
import logging
from ....utils.exception import ZeroFlowError
from scipy.constants import atm, R, zero_Celsius

MOLAR_MASS_AIR = 28.97  # g/mol

FLOW_EQUATION_CONSTANT = math.pi * math.sqrt(R/16./MOLAR_MASS_AIR*1000)


def calculate_pipeline_average_temperature(t_ambient, t1, t2):
    """

    :param t_ambient: Ambient temperature
    :param t1: Inlet temperature
    :param t2: Outlet temperature
    :return: Pipeline average temperature
    """
    if t1 == t_ambient or t2 == t_ambient or t1 == t2:
        return t1
    else:
        return t_ambient + (t1 - t2) / math.log((t1 - t_ambient) / (t2 - t_ambient))


def calculate_pipeline_average_pressure(p1, p2):
    """

    :param p1:
    :param p2:
    :return:
    """
    return 2./3. * ((p1 + p2) - (p1 * p2) / (p1 + p2))


def calculate_height_difference_correction(h1, h2, z, p, t, sg):
    """
    Calculate the effect on the flow rate calculation caused by the height difference between pipeline inlet and outlet
    :param h1: Inlet height [m]
    :param h2: Outlet height [m]
    :param z: Gas mixture compressibility factor [-]
    :param p: Pipeline average pressure [Pa]
    :param t: Pipeline average temperature [K]
    :param sg: Gas mixture specific gravity [-]
    :return:
    """
    if h1 == h2:
        return 0
    else:
        return 0.06835 * sg * (h2 - h1) * p ** 2 / (z * t)

    # try:
    #     return 0.06835 * sg * (h2 - h1) * p ** 2 / (z * t)
    # except:
    #     print("Error calculating the slope correction!")


def calculate_pipeline_physical_characteristic(d, length, t_avg, f, eta):
    tb = zero_Celsius + 15  # Temperature base, 15 Celsius
    pb = 1 * atm  # Pressure base, 1 atm

    return (FLOW_EQUATION_CONSTANT * tb / pb) * (d ** 2.5) * \
           ((1 / (length * t_avg * f)) ** 0.5) * eta


def determine_flow_direction(p1, p2, h_diff_correction):
    tmp = calculate_energy_difference(p1, p2, h_diff_correction)
    if tmp > 0:
        return 1
    elif tmp < 0:
        return -1
    else:
        raise ZeroFlowError('Got zero flow in pipeline.')


def calculate_volumetric_flow_rate(p1, p2, h_diff_correction, pipe_physical_char, z, sg):
    f_direction = determine_flow_direction(p1, p2, h_diff_correction)
    return f_direction * (abs(p1 ** 2 - p2 ** 2 - h_diff_correction)/z/sg) ** (1 / 2) * pipe_physical_char


def calculate_mass_flow_rate(q, rho_standard):
    """

    :param q: Volumetric flow rate at the standard condition [sm3/s]
    :param rho_standard: Gas mixture density at the standard condition [kg/sm3]
    :return: Mass flow rate [kg/s]
    """
    return rho_standard * q


def calculate_energy_difference(p1, p2, h_diff_correction):
    return p1 ** 2 - p2 ** 2 - h_diff_correction


def volumetric_flow_rate_first_order_derivative(pipe_physical_char, energy_difference, p):
    return pipe_physical_char * (abs(energy_difference)) ** (-0.5)


def mass_flow_rate_first_order_derivative(rho_standard, q_first_order_derivative):
    return rho_standard * q_first_order_derivative
