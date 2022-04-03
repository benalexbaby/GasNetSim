#  ******************************************************************************
#   Copyright (c) 2020.
#   Developed by Yifei Lu
#  ******************************************************************************

# from GasNetSim.components.gas_mixture.thermo.thermo import Mixture
import math
from collections import OrderedDict
import matplotlib.pyplot as plt
import os
from pathlib import Path
from thermo import Mixture


def calc_papay_compressibility(pressure, temperature, h2_fraction):
    # pressure = pressure / 101.325
    critical_pressure_ch4 = 45.79  # bar
    critical_temperature_ch4 = 190.8  # K
    critical_pressure_h2 = 12.8  # bar
    critical_temperature_h2 = 33.2  # K

    h2_papay_z = 1 - 3.52 * (pressure / critical_pressure_h2) * math.exp(-2.26 * (temperature / critical_temperature_h2)) + \
           0.274 * (pressure/critical_pressure_h2)**2 * math.exp(-1.878 * (temperature / critical_temperature_h2))

    ch4_papay_z = 1 - 3.52 * (pressure / critical_pressure_ch4) * math.exp(-2.26 * (temperature / critical_temperature_ch4)) + \
           0.274 * (pressure/critical_pressure_ch4)**2 * math.exp(-1.878 * (temperature / critical_temperature_ch4))

    return h2_papay_z * h2_fraction + ch4_papay_z * (1-h2_fraction)

def calc_aga_compressibility(pressure, temperature, h2_fraction):
    # pressure = pressure / 101.325
    critical_pressure_ch4 = 45.79  # bar
    critical_temperature_ch4 = 190.8  # K
    critical_pressure_h2 = 12.8  # bar
    critical_temperature_h2 = 33.2  # K

    h2_aga_z = 1 + 0.257 * (pressure / critical_pressure_h2) - \
                 0.533 * (pressure/critical_pressure_h2)/(temperature / critical_temperature_h2)

    ch4_aga_z = 1 + 0.257 * (pressure / critical_pressure_ch4) - \
                 0.533 * (pressure/critical_pressure_ch4)/(temperature / critical_temperature_ch4)

    return h2_aga_z * h2_fraction + ch4_aga_z * (1-h2_fraction)


save_folder = Path("./figures/")

gas_comp = OrderedDict([('methane', 1.0),
                        ('hydrogen', 0.0)])

gas_mix_z = list()
gas_mix_sg = list()
gas_mix_papay_z = list()
# gas_mix_aga_z = list()

# gas_mix_gerg_z = list()
#
# with open(Path('../../gas_mixture/GERG2008/output.txt'), 'r') as f:
#     for line in f:
#         gas_mix_gerg_z.append(float(line))

while gas_comp['methane'] > -0.01:
    gas_mixture = Mixture(P=20 * 101325, T=300, zs=gas_comp)
    print("The gas mixture compressibility is {:f}.".format(gas_mixture.Z))
    print("The gas relative density is {:f}.".format(gas_mixture.SG))
    gas_mix_z.append(gas_mixture.Z)
    gas_mix_sg.append(gas_mixture.SG)
    gas_mix_papay_z.append(calc_papay_compressibility(20, 300, gas_comp['hydrogen']))
    # gas_mix_aga_z.append(calc_aga_compressibility(50, 300, gas_comp['hydrogen']))
    gas_comp['methane'] -= 0.01
    gas_comp['hydrogen'] += 0.01

h2_z = list()
for x in range(1, 100):
    hydrogen_compressibility = Mixture(zs={'hydrogen':1.0}, P=x * 101325, T=300).Z
    h2_z.append(hydrogen_compressibility)

plt.style.use("ieeetrans")
# import matplotlib
#
# matplotlib.rcParams['mathtext.fontset'] = 'custom'
# matplotlib.rcParams['mathtext.rm'] = 'Times New Roman'
# matplotlib.rcParams['mathtext.it'] = 'Times New Roman:italic'
# matplotlib.rcParams['mathtext.bf'] = 'Times New Roman:bold'

plt.figure()
plt.plot(gas_mix_z, label="PREOS")
# plt.plot([0, 100], [gas_mix_z[0], gas_mix_z[100]], label="linear model")
plt.plot(gas_mix_papay_z, label="linear Papay's equation")
# plt.plot(gas_mix_gerg_z, label="GERG-2008")
# plt.plot(gas_mix_aga_z, label="linear AGA")
# plt.title("Compressibility considering hydrogen injection")
plt.xlabel("Hydrogen concentration [%]")
plt.ylabel("Gas mixture compressibility")
plt.legend()
# plt.savefig(save_folder.joinpath("Z_h2_injection.pdf"))
plt.show()

plt.figure()
plt.plot(gas_mix_sg, label="PREOS")
plt.plot([0, 100], [gas_mix_sg[0], gas_mix_sg[100]], label="linear model")
# plt.title("Specific gravity considering hydrogen injection")
plt.xlabel("Hydrogen concentration [%]")
plt.ylabel("Gas mixture specific gravity")
plt.legend()
# plt.savefig(save_folder.joinpath("SG_h2_injection.pdf"))
plt.show()

from thermo import Mixture

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

gas_mixture = Mixture(zs=gas_comp, T=288.15, P=50 * 101325)
