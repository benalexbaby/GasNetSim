from numpy.testing import assert_almost_equal
from collections import OrderedDict
from scipy.constants import bar
import pandas as pd
import math
import matplotlib.pyplot as plt

from GasNetSim.components.utils.gas_mixture.GERG2008.gerg2008 import *
import cantera as ct

def heating_value(fuel):
    """Returns the LHV and HHV for the specified fuel"""
    gas.TP = 298, ct.one_atm
    gas.set_equivalence_ratio(1.0, fuel, "O2:1.0")
    h1 = gas.enthalpy_mass
    Y_fuel = gas[fuel].Y[0]

    # complete combustion products
    X_products = {
        "CO2": gas.elemental_mole_fraction("C"),
        "H2O": 0.5 * gas.elemental_mole_fraction("H"),
        "N2": 0.5 * gas.elemental_mole_fraction("N"),
    }

    gas.TPX = None, None, X_products
    Y_H2O = gas["H2O"].Y[0]
    h2 = gas.enthalpy_mass
    LHV = -(h2 - h1) / Y_fuel / 1e6
    HHV = -(h2 - h1 + (h_liquid - h_gas) * Y_H2O) / Y_fuel / 1e6
    return LHV, HHV

if __name__ == '__main__':
    fuels = ["CH4", "C2H6", "C3H8", "H2", "CO"]
    molar_density = [0.6484926588314163, 1.2228248480248802, 1.8085351993461924, 0.04065274708618872 * 2, 1.1308810162414191]
    gas = ct.Solution("gri30.yaml")
    water = ct.Water()
    # Set liquid water state, with vapor fraction x = 0
    water.TQ = 298, 0
    h_liquid = water.h
    # Set gaseous water state, with vapor fraction x = 1
    water.TQ = 298, 1
    h_gas = water.h
    # cantera MJ/kg - MJ/m3
    LHV_cantera_mass = []
    HHV_cantera_mass = []
    LHV_cantera_vol = []
    HHV_cantera_vol = []
    for i in range(len(fuels)):
        LHV, HHV = heating_value(fuels[i])
        LHV_cantera_mass.append(LHV)
        HHV_cantera_mass.append(HHV)
        LHV_cantera_vol.append(LHV * molar_density[i])
        HHV_cantera_vol.append(HHV * molar_density[i])
    # gerg2008 kJ/kg - kJ/m3
    LHV_gerg2008_mass = []
    HHV_gerg2008_mass = []
    LHV_gerg2008_vol = []
    HHV_gerg2008_vol = []
    gas_comp = {'methane': 1, 'ethane': 1, 'propane': 1, 'hydrogen': 1, 'carbon monoxide': 1}
    for key, value in gas_comp.items():
        gas_mixture = GasMixtureGERG2008(P=1 * bar, T=298, composition={key: value})
        HHV = gas_mixture.CalculateHeatingValue({key: value}, hhv=True, parameter='mass')
        LHV = gas_mixture.CalculateHeatingValue({key: value}, hhv=False, parameter='mass')
        LHV_gerg2008_mass.append(LHV)
        HHV_gerg2008_mass.append(HHV)
        HHV = gas_mixture.CalculateHeatingValue({key: value}, hhv=True, parameter='vol')
        LHV = gas_mixture.CalculateHeatingValue({key: value}, hhv=False, parameter='vol')
        LHV_gerg2008_vol.append(LHV)
        HHV_gerg2008_vol.append(HHV)
    LHV_gerg2008_mass = [x / 1000 for x in LHV_gerg2008_mass]
    LHV_gerg2008_vol = [x / 1000 for x in LHV_gerg2008_vol]
    HHV_gerg2008_mass = [x / 1000 for x in HHV_gerg2008_mass]
    HHV_gerg2008_vol = [x / 1000 for x in HHV_gerg2008_vol]

    LHV_df_mass = pd.DataFrame(
        {'molecule': fuels,
         'gerg2008 [LHV MJ/kg]': LHV_gerg2008_mass,
         'cantera [LHV MJ/kg]': LHV_cantera_mass,
         })

    print(LHV_df_mass)

    HHV_df_mass = pd.DataFrame(
        {'molecule': fuels,
         'gerg2008 [HHV MJ/kg]': HHV_gerg2008_mass,
         'cantera [HHV MJ/kg]': HHV_cantera_mass
         })

    print(HHV_df_mass)

    LHV_df_vol = pd.DataFrame(
        {'molecule': fuels,
         'gerg2008 [LHV MJ/m3]': LHV_gerg2008_vol,
         'cantera [LHV MJ/m3]': LHV_cantera_vol,
         })

    print(LHV_df_vol)

    HHV_df_vol = pd.DataFrame(
        {'molecule': fuels,
         'gerg2008 [HHV MJ/m3]': HHV_gerg2008_vol,
         'cantera [HHV MJ/m3]': HHV_cantera_vol,
         })

    print(HHV_df_vol)

    # assert almost equal
    assert_almost_equal(LHV_cantera_mass, LHV_gerg2008_mass, decimal=2)
    assert_almost_equal(HHV_cantera_mass, HHV_gerg2008_mass, decimal=2)
    assert_almost_equal(LHV_cantera_vol, LHV_gerg2008_vol, decimal=2)
    assert_almost_equal(HHV_cantera_vol, HHV_gerg2008_vol, decimal=2)

    cantera_list = LHV_cantera_mass + HHV_cantera_mass + LHV_cantera_vol + HHV_cantera_vol
    gerg2008_list = LHV_gerg2008_mass + HHV_gerg2008_mass + LHV_gerg2008_vol + HHV_gerg2008_vol
    assert_almost_equal(cantera_list, gerg2008_list, decimal=2)

    # print("fuel   LHV (MJ/kg)   HHV (MJ/kg)")
    # for fuel in fuels:
    #
    #     print(f"{fuel:8s} {LHV:7.7f}      {HHV:7.7f}")