#  #!/usr/bin/env python
#  -*- coding: utf-8 -*-
#  ******************************************************************************
#    Copyright (c) 2022.
#    Developed by Yifei Lu
#    Last change on 3/9/22, 11:18 AM
#    Last change by yifei
#   *****************************************************************************
import matplotlib.pyplot as plt
import numpy as np
import math


def calc_beta(ul, qm, cp, d):
    """

    :param ul: Heat transfer coefficient [W m^(-2) K^(-1)]
    :param qm: Mass flow [kg/s]
    :param cp: Specific heat capacity at constant pressure [J/(kg K)]
    :param d: pipeline diameter [m]
    :return:
    """
    return ul/(qm*cp)*np.pi*d


def calc_gamma(mu_jt, z, R, f, qm, p_a, D):
    g = 9.81
    A = math.pi * (D/2)**2
    return mu_jt * z * R * f * qm * abs(qm) / (2*p_a*D*A**2)


def outlet_temp(beta, gamma, Ts, L, T1):
    return beta/(beta+gamma)*(Ts - Ts * math.exp(-(beta+gamma)*L)) + T1*math.exp(-(beta+gamma)*L)
