#   #!/usr/bin/env python
#   -*- coding: utf-8 -*-
#   ******************************************************************************
#     Copyright (c) 2022.
#     Developed by Yifei Lu
#     Last change on 3/22/22, 12:31 AM
#     Last change by yifei
#    *****************************************************************************
import matplotlib.pyplot as plt
import GasNetSim as Gns
from GasNetSim.components.utils.gas_mixture.heating_value import *

if __name__ == "__main__":
    from collections import OrderedDict

    gas_comp = OrderedDict([('methane', 1.0),
                            ('hydrogen', 0.0)])
    step_size = 0.2

    plt.figure()
    for temperature in range(240, 340, 20):
        gas_mixture = Gns.Mixture(zs=gas_comp, T=temperature, P=50 * 101325)
        HHV = calc_heating_value(gas_mixture, heating_value_type='HHV')
        mass_flow_rate = 1000 * 1e6 / HHV
        temp = list()
        beta = calc_beta(3.69, mass_flow_rate, gas_mixture.Cp, d=0.5)
        for l in range(1, 150):
            temp.append(outlet_temp(beta=beta,
                                    gamma=calc_gamma(gas_mixture.JT, gas_mixture.Z, gas_mixture.R_specific,
                                                     f=0.01, qm=mass_flow_rate, p_a=50 * 101325, D=0.5),
                                    Ts=288.15, L=l * 1000, T1=temperature))
        plt.plot(temp, label=f'T1 = {temperature} K')
        plt.legend(title='Inlet temperature [K]')
        plt.xlabel('Pipe length [km]')
        plt.ylabel('Temperature [K]')
    # plt.savefig('figures/temperature_different_t1.png')
    plt.show()

    plt.figure()
    while gas_comp['methane'] >= -0.01:
        temp = list()
        gas_mixture = Mixture(zs=gas_comp, T=300, P=50 * 101325)
        HHV = calc_heating_value(gas_mixture)
        mass_flow_rate = 1000*1e6/HHV
        beta = calc_beta(3.69, mass_flow_rate, gas_mixture.Cp, d=0.5)
        for l in range(1, 150):
            temp.append(outlet_temp(beta=beta,
                                    gamma=calc_gamma(gas_mixture.JT, gas_mixture.Z, gas_mixture.R_specific,
                                                    f=0.01, qm=mass_flow_rate, p_a=50*101325, D=0.5),
                                    Ts=288.15, L=l*1000, T1=300))
        plt.plot(temp, label=f'{int(gas_comp["hydrogen"]*100)}%')
        plt.legend(title='Hydrogen concentration')
        plt.xlabel('Pipe length [km]')
        plt.ylabel('Temperature [K]')
        gas_comp['methane'] -= step_size
        gas_comp['hydrogen'] += step_size
    # plt.savefig('figures/temperature_different_composition.png')
    plt.show()
