#  #!/usr/bin/env python
#  -*- coding: utf-8 -*-
#  ******************************************************************************
#    Copyright (c) 2021.
#    Developed by Yifei Lu
#    Last change on 12/21/21, 4:58 PM
#    Last change by yifei
#   *****************************************************************************
from collections import OrderedDict
import cantera as ct

# from .thermo.thermo import Mixture
from thermo import Mixture
from .gas_mixture import *

GAS = ct.Solution('gri30.yaml')


def get_mole_fraction(gas_mixture):
    """
    Get mole fraction of the gas composition at node
    :return: Gas mole fraction
    """
    thermo_mixture = Mixture(zs=gas_mixture.composition,
                             T=gas_mixture.temperature,
                             P=gas_mixture.pressure)
    mole_fraction = dict()
    for i in range(len(thermo_mixture.zs)):
        gas = thermo_mixture.formulas[i]
        mole_fraction[gas] = thermo_mixture.zs[i]
    return mole_fraction


def calc_heating_value(gas_mixture, heating_value_type='HHV'):
    """ Returns the LHV and HHV for the specified fuel """
    GAS.TPX = gas_mixture.temperature, gas_mixture.pressure, get_mole_fraction(gas_mixture)
    GAS.set_equivalence_ratio(1.0, get_mole_fraction(gas_mixture), 'O2:1.0')
    h1 = GAS.enthalpy_mass
    Y_fuel = 1 - GAS['O2'].Y[0]

    # complete combustion products
    Y_products = {'CO2': GAS.elemental_mole_fraction('C'),
                  'H2O': 0.5 * GAS.elemental_mole_fraction('H'),
                  'N2': 0.5 * GAS.elemental_mole_fraction('N')}

    GAS.TPX = None, None, Y_products
    Y_H2O = GAS['H2O'].Y[0]
    h2 = GAS.enthalpy_mass
    # Calculate LHV
    LHV = -(h2-h1)/Y_fuel
    # Calculate HHV
    water = ct.Water()
    water.TQ = 298, 0  # Set liquid water state, with vapor fraction x = 0
    h_liquid = water.h
    water.TQ = 298, 1  # Set gaseous water state, with vapor fraction x = 1
    h_gas = water.h
    HHV = -(h2-h1 + (h_liquid-h_gas) * Y_H2O)/Y_fuel
    if heating_value_type == 'HHV':
        return HHV
    else:
        return LHV


if __name__ == "__main__":
    gas_comp = OrderedDict([('methane', 1.0),
                            ('hydrogen', 0.0)])

    print('fuel   LHV (MJ/kg)   HHV (MJ/kg)')
    while gas_comp['methane'] > -0.01:
        gas_mixture = GasMixture(P=70 * 101325, T=300, zs=gas_comp)
        standard_density = Mixture(P=101325, T=288.15, zs=gas_comp).rho
        gas_comp['methane'] -= 0.01
        gas_comp['hydrogen'] += 0.01
        HHV = calc_heating_value(gas_mixture)
        LHV = calc_heating_value(gas_mixture, heating_value_type='LHV')
        print(' {:7.3f}      {:7.3f}'.format(LHV/1e6, HHV/1e6))
