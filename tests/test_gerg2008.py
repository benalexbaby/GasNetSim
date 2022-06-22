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


# NATURAL_GAS = [0, 0.94, 0.05, 0.01] + [0] * 18
# gas = GasMixtureGERG2008(P=101325, T=288.15, x=NATURAL_GAS)
P = 20 * 101325
T = 300

while NATURAL_GAS['methane']> 0:
    gas = GasMixtureGERG2008(P=P, T=T, composition=NATURAL_GAS)
    print(gas.Z)
    NATURAL_GAS['methane'] -= 0.01
    NATURAL_GAS['hydrogen'] += 0.01
