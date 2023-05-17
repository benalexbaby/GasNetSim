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
from collections import OrderedDict
from scipy.constants import bar
import pandas as pd
import math

# test over natural gas
NATURAL_GAS = OrderedDict([('methane', 0.947),
                           ('ethane', 0.042),
                           ('propane', 0.002),
                           ('isobutane', 0.0002),
                           ('butane', 0.0002),
                           ('isopentane', 0.0001),
                           ('pentane', 0.0001),
                           ('hexane', 0.0001),
                           ('nitrogen', 0.005),
                           ('carbon dioxide', 0.003),
                           ('oxygen', 0.0001),
                           ('hydrogen', 0.0002)])

gas_mixture = GasMixtureGERG2008(P_Pa=1 * bar, T_K=298, composition=NATURAL_GAS)
print("Molar Mass [g/mol] = "+str(gas_mixture.MolarMass))
print("Z [-] = "+str(gas_mixture.Z))
print("Isochoric Heat Capacity [J/mol-K] = "+str(gas_mixture.Cv))
print("Isobaric Heat Capacity [J/mol-K] = "+str(gas_mixture.Cp))
print("Joule-Thomson coefficient [K/kPa] = "+str(gas_mixture.JT))

# test over methane and hydrogen blending
z_list = []
h_list = []

gas_comp = OrderedDict([('methane', 1.0), ('hydrogen', 0.0)])
h_list.append(gas_comp['hydrogen'])
gas_mixture = GasMixtureGERG2008(P_Pa=1 * bar, T_K=298, composition=gas_comp)

Z = gas_mixture.Z
z_list.append(Z)

print('hydrogen content [%]      Z [-] ')
print('{:7.3f}      {:7.3f}'.format(gas_comp['hydrogen']*100,Z))

while gas_comp['methane'] >= 0:

    gas_comp['methane'] -= 0.01
    gas_comp['hydrogen'] += 0.01
    h_list.append(gas_comp['hydrogen'])

    gas_mixture = GasMixtureGERG2008(P_Pa=1 * bar, T_K=298, composition=gas_comp)

    Z = gas_mixture.Z
    z_list.append(Z)

    print('{:7.3f}      {:7.3f}'.format(gas_comp['hydrogen']*100,Z))
