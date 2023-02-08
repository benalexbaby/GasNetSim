#   #!/usr/bin/env python
#   -*- coding: utf-8 -*-
#   ******************************************************************************
#     Copyright (c) 2023.
#     Developed by Yifei Lu
#     Last change on 1/6/23, 3:03 PM
#     Last change by yifei
#    *****************************************************************************
from GasNetSim.components.utils.gas_mixture.GERG2008.gerg2008 import *
from scipy.constants import bar
from numpy.testing import assert_almost_equal


def test_gerg_function_with_original_cpp_implementation():
    nist_gas_mixture = {}
    a = ['methane', 'nitrogen', 'carbon dioxide', 'ethane', 'propane', 'isobutane',
     'butane', 'isopentane', 'pentane', 'hexane', 'heptane', 'octane', 'nonane',
     'decane', 'hydrogen', 'oxygen', 'carbon monoxide', 'water', 'hydrogen sulfide',
     'helium', 'argon']
    b = [0.77824, 0.02, 0.06, 0.08, 0.03, 0.0015, 0.003, 0.0005, 0.00165, 0.00215, 0.00088, 0.00024, 0.00015, 0.00009,
         0.004, 0.005, 0.002, 0.0001, 0.0025, 0.007, 0.001]
    for i in range(21):
        nist_gas_mixture[a[i]] = b[i]

    gas_mixture = GasMixtureGERG2008(500 * bar, 400, nist_gas_mixture)

    nist_results = {"Molar mass [g/mol]": 20.5427445016,
                    "Molar density [mol/l]": 12.79828626082062,
                    "Pressure [kPa]": 50000.0,
                    "Compressibility factor": 1.174690666383717,
                    "d(P)/d(rho) [kPa/(mol/l)]": 7000.694030193327,
                    "d^2(P)/d(rho)^2 [kPa/(mol/l)^2]": 1129.526655214841,
                    "d(P)/d(T) [kPa/K]": 235.9832292593096,
                    "Energy [J/mol]": -2746.492901212530,
                    "Enthalpy [J/mol]": 1160.280160510973,
                    "Entropy [J/mol-K]": -38.57590392409089,
                    "Isochoric heat capacity [J/mol-K]": 39.02948218156372,
                    "Isobaric heat capacity [J/mol-K]": 58.45522051000366,
                    "Speed of sound [m/s]": 714.4248840596024,
                    "Gibbs energy [J/mol]": 16590.64173014733,
                    "Joule-Thomson coefficient [K/kPa]": 7.155629581480913E-05,
                    "Isentropic exponent": 2.683820255058032}

    mm = gas_mixture.MolarMass
    D = gas_mixture.MolarDensity
    Z = gas_mixture.Z
    dPdD = gas_mixture.dPdD
    d2PdD2 = gas_mixture.d2PdD2
    dPdT = gas_mixture.dPdT
    U = gas_mixture.energy
    H = gas_mixture.enthalpy
    S = gas_mixture.entropy
    Cv = gas_mixture.Cv
    Cp = gas_mixture.Cp
    W = gas_mixture.c
    G = gas_mixture.gibbs_energy
    JT = gas_mixture.JT
    Kappa = gas_mixture.isentropic_exponent

    my_results = {"Molar mass [g/mol]": mm,
                  "Molar density [mol/l]": D,
                  "Pressure [kPa]": gas_mixture.P,
                  "Compressibility factor": Z,
                  "d(P)/d(rho) [kPa/(mol/l)]": dPdD,
                  "d^2(P)/d(rho)^2 [kPa/(mol/l)^2]": d2PdD2,
                  "d(P)/d(T) [kPa/K]": dPdT,
                  "Energy [J/mol]": U,
                  "Enthalpy [J/mol]": H,
                  "Entropy [J/mol-K]": S,
                  "Isochoric heat capacity [J/mol-K]": Cv,
                  "Isobaric heat capacity [J/mol-K]": Cp,
                  "Speed of sound [m/s]": W,
                  "Gibbs energy [J/mol]": G,
                  "Joule-Thomson coefficient [K/kPa]": JT*1e3,
                  "Isentropic exponent": Kappa}

    assert_almost_equal(nist_results["Molar mass [g/mol]"], my_results["Molar mass [g/mol]"])
    assert_almost_equal(nist_results["Molar density [mol/l]"], my_results["Molar density [mol/l]"])
    assert_almost_equal(nist_results["Pressure [kPa]"], my_results["Pressure [kPa]"])
    assert_almost_equal(nist_results["Compressibility factor"], my_results["Compressibility factor"])
    assert_almost_equal(nist_results["d(P)/d(rho) [kPa/(mol/l)]"], my_results["d(P)/d(rho) [kPa/(mol/l)]"])
    assert_almost_equal(nist_results["d^2(P)/d(rho)^2 [kPa/(mol/l)^2]"], my_results["d^2(P)/d(rho)^2 [kPa/(mol/l)^2]"])
    assert_almost_equal(nist_results["d(P)/d(T) [kPa/K]"], my_results["d(P)/d(T) [kPa/K]"])
    assert_almost_equal(nist_results["Energy [J/mol]"], my_results["Energy [J/mol]"])
    assert_almost_equal(nist_results["Enthalpy [J/mol]"], my_results["Enthalpy [J/mol]"])
    assert_almost_equal(nist_results["Entropy [J/mol-K]"], my_results["Entropy [J/mol-K]"])
    assert_almost_equal(nist_results["Isochoric heat capacity [J/mol-K]"],
                        my_results["Isochoric heat capacity [J/mol-K]"])
    assert_almost_equal(nist_results["Isobaric heat capacity [J/mol-K]"],
                        my_results["Isobaric heat capacity [J/mol-K]"])
    assert_almost_equal(nist_results["Speed of sound [m/s]"], my_results["Speed of sound [m/s]"])
    assert_almost_equal(nist_results["Gibbs energy [J/mol]"], my_results["Gibbs energy [J/mol]"])
    assert_almost_equal(nist_results["Joule-Thomson coefficient [K/kPa]"],
                        my_results["Joule-Thomson coefficient [K/kPa]"])
    assert_almost_equal(nist_results["Isentropic exponent"], my_results["Isentropic exponent"])
