#   #!/usr/bin/env python
#   -*- coding: utf-8 -*-
#   ******************************************************************************
#     Copyright (c) 2022.
#     Developed by Yifei Lu
#     Last change on 4/3/22, 1:50 PM
#     Last change by yifei
#    *****************************************************************************
from GasNetSim.components.utils.gas_mixture.typical_mixture_composition import NATURAL_GAS
from GasNetSim.components.utils.gas_mixture.GERG2008 import *

from scipy.constants import bar


def test_gerg2008_gas_mixture_properties():
    P = 20 * bar
    T = 300

    while NATURAL_GAS['methane'] > 0:
        gas = GasMixtureGERG2008(P_Pa=P, T_K=T, composition=NATURAL_GAS)
        NATURAL_GAS['methane'] -= 0.01
        NATURAL_GAS['hydrogen'] += 0.01

