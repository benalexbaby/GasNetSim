# """
# Version 2.01 of routines for the calculation of thermodynamic
# properties from the AGA 8 Part 2 GERG-2008 equation of state.
# April, 2017
#
# Written by Eric W. Lemmon
# Applied Chemicals and Materials Division
# National Institute of Standards and Technology (NIST)
# Boulder, Colorado, USA
# Eric.Lemmon@nist.gov
# 303-497-7939
#
# Python translation by Yifei Lu
#
# IEK-10: Energy Systems Engineering
# Institute of Energy and Climate Research, Forschungszentrum Juelich
# Juelich, NRW, Germany
# yifei.lu@fz-juelich.de
#
# Other contributors:
# Volker Heinemann, RMG Messtechnik GmbH
# Jason Lu, Thermo Fisher Scientific
# Ian Bell, NIST
# """
#
# """
# The publication for the AGA 8 equation of state is available from AGA
#   and the Transmission Measurement Committee.
#
# The GERG-2008 equation of state was developed by Oliver Kunz and Wolfgang Wagner;
#
# Kunz, O. and Wagner, W.
# The GERG-2008 Wide-Range Equation of State for Natural Gases and Other Mixtures;
# An Expansion of GERG-2004
# J. Chem. Eng. Data, 57(11):3032-3091, 2012.
#
# Kunz, O., Klimeck, R., Wagner, W., and Jaeschke, M.
# The GERG-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
# GERG Technical Monograph 15
# Fortschr.-Ber. VDI, Reihe 6, Nr. 557, VDI Verlag, Düsseldorf, 2007.
# http://www.gerg.eu/public/uploads/files/publications/technical_monographs/tm15_04.pdf
# """
#
# """
# Subroutines contained here for property calculations:
# ***** Subroutine SetupGERG must be called once before calling other routines. ******
# Sub MolarMassGERG(x, Mm)
# Sub PressureGERG(T, D, x, P, Z)
# Sub DensityGERG(iFlag, T, P, x, D, ierr, herr)
# Sub PropertiesGERG(T, D, x, P, Z, dPdD, d2PdD2, d2PdTD, dPdT, U, H, S, Cv, Cp, W, G, JT, Kappa)
# Sub SetupGERG()
# """
#
# """
# The compositions in the x() array use the following order and must be sent as mole fractions:
#     1 - Methane
#     2 - Nitrogen
#     3 - Carbon dioxide
#     4 - Ethane
#     5 - Propane
#     6 - Isobutane
#     7 - n-Butane
#     8 - Isopentane
#     9 - n-Pentane
#    10 - n-Hexane
#    11 - n-Heptane
#    12 - n-Octane
#    13 - n-Nonane
#    14 - n-Decane
#    15 - Hydrogen
#    16 - Oxygen
#    17 - Carbon monoxide
#    18 - Water
#    19 - Hydrogen sulfide
#    20 - Helium
#    21 - Argon
#
# For example, a mixture (in moles) of 94% methane, 5% CO2, and 1% helium would be (in mole fractions):
# x(1)=0.94, x(3)=0.05, x(20)=0.01
# """
#
# """
# // Function prototypes (not exported)
# static void Alpha0GERG(const double T, const double D, const std::vector<double> &x, double a0[3]);
# static void AlpharGERG(const int itau, const int idelta, const double T, const double D, const std::vector<double> &x, double ar[4][4]);
# static void PseudoCriticalPointGERG(const std::vector<double> &x, double &Tcx, double &Dcx);
# static void tTermsGERG(const double lntau, const std::vector<double> &x);
# """
from collections import OrderedDict

#   #!/usr/bin/env python
#   -*- coding: utf-8 -*-
#   ******************************************************************************
#     Copyright (c) 2022.
#     Developed by Yifei Lu
#     Last change on 4/3/22, 1:38 PM
#     Last change by yifei
#    *****************************************************************************

import numpy as np
import math
from collections import Counter
from .setup import *
from copy import deepcopy

from GasNetSim.components.utils.gas_mixture.GERG2008.setup import *


def Tanh(xx):
    return (math.exp(xx) - math.exp(-xx)) / (math.exp(xx) + math.exp(-xx))


def Sinh(xx):
    return (math.exp(xx) - math.exp(-xx)) / 2


def Cosh(xx):
    return (math.exp(xx) + math.exp(-xx)) / 2


"""
The compositions in the x() array use the following order and must be sent as mole fractions:
    1 - Methane
    2 - Nitrogen
    3 - Carbon dioxide
    4 - Ethane
    5 - Propane
    6 - Isobutane
    7 - n-Butane
    8 - Isopentane
    9 - n-Pentane
   10 - n-Hexane
   11 - n-Heptane
   12 - n-Octane
   13 - n-Nonane
   14 - n-Decane
   15 - Hydrogen
   16 - Oxygen
   17 - Carbon monoxide
   18 - Water
   19 - Hydrogen sulfide
   20 - Helium
   21 - Argon
"""


class GasMixtureGERG2008:
    def __init__(self, P_Pa, T_K, composition):
        # Input parameters
        self.dPdT = None
        self.d2PdD2 = None
        self.dPdD = None
        self.P = P_Pa / 1000  # Pa -> kPa

        self.T = T_K
        self.x = self.CovertCompositionGERG(composition=composition)  # gas composition

        # Calculated properties
        self.MolarMass = self.MolarMassGERG()
        # self.MolarDensity = self.DensityGERG(iFlag=0)[2]
        self.MolarDensity = 0
        self.rho = 0
        self.SG = 1
        self.Z = 1
        self.energy = 0
        self.enthalpy = 0
        self.entropy = 0
        self.Cv_molar = 0  # molar isochoric heat capacity [J/mol-K]
        self.Cp_molar = 0  # molar isobaric heat capacity [J/mol-K]
        self.Cp = 0  # isochoric heat capacity [J/kg-K]
        self.Cv = 0  # isobaric heat capacity [J/kg-K]
        self.c = 0  # speed of sound [m/s]
        self.gibbs_energy = 0  # Gibbs energy [J/mol]

        self.JT = 1  # Joule-Thomson coefficient [K/Pa]
        self.isentropic_exponent = 0  # Isentropic exponent

        self.R_specific = 0
        self.viscosity = 2e-4  # TODO add function

        self.PropertiesGERG()

        self.HHV = self.CalculateHeatingValue(comp=composition, hhv=True, parameter="volume")

    def CalculateHeatingValue(self, comp, hhv, parameter):

        dict_components = {'methane': {'C': 1, 'H': 4},
                           'nitrogen': {'N': 2},
                           'carbon dioxide': {'C': 1, 'O': 2},
                           'ethane': {'C': 2, 'H': 6},
                           'propane': {'C': 3, 'H': 8},
                           'isobutane': {'C': 4, 'H': 10},
                           'n-butane': {'C': 4, 'H': 10},
                           'isopentane': {'C': 5, 'H': 12},
                           'n-pentane': {'C': 5, 'H': 12},
                           'n-hexane': {'C': 6, 'H': 14},
                           'n-heptane': {'C': 7, 'H': 16},
                           'n-octane': {'C': 8, 'H': 18},
                           'n-nonane': {'C': 9, 'H': 20},
                           'n-decane': {'C': 10, 'H': 22},
                           'hydrogen': {'H': 2},
                           'oxygen': {'O': 2},
                           'carbon monoxide': {'C': 1, 'O': 1},
                           'water': {'H': 2, 'O': 1},
                           'hydrogen sulfide': {'H': 2, 'S': 1},
                           'helium': {'He': 1},
                           'argon': {'Ar': 1}}
                           # 'SO2': {'S': 1, 'O': 2}}
        # 298 K
        # dict_enthalpy_mole = {'methane': -74602.416533355,
        #                       'nitrogen': 0.0,
        #                       'carbon dioxide': -393517.79827154,
        #                       'ethane': -83856.2627150042,
        #                       'propane': -103861.117481869,
        #                       'isobutane': -135360.0,
        #                       'n-butane': -125849.99999999999,
        #                       'isopentane': -178400.0,
        #                       'n-pentane': -173500.0,
        #                       'n-hexane': -198490.0,
        #                       'n-heptane': -223910.0,
        #                       'n-hctane': -249730.0,
        #                       'n-nonane': -274700.0,
        #                       'n-decane': -300900.0,
        #                       'hydrogen': 0.0,
        #                       'oxygen': -4.40676212751828,
        #                       'carbon monoxide': -110525.0,
        #                       'water': -241833.418361837,
        #                       'hydrogen sulfide': -20600.0,
        #                       'helium': 0.0,
        #                       'argon': 0.0,
        #                       'carbon': 0.0}
        #                       # 'H': 218000.0,
        #                       # 'O': 249190.0,
        #                       # 'SO2': -296840.0}

        # 273 K
        dict_enthalpy_mole = {'methane': -75483.51423273719,
                              'nitrogen': 0.0,
                              'carbon dioxide': -394431.82606764464,
                              'ethane': -83856.2627150042,
                              'propane': -103861.117481869,
                              'isobutane': -135360.0,
                              'n-butane': -125849.99999999999,
                              'isopentane': -178400.0,
                              'n-pentane': -173500.0,
                              'n-hexane': -198490.0,
                              'n-heptane': -223910.0,
                              'n-hctane': -249730.0,
                              'n-nonane': -274700.0,
                              'n-decane': -300900.0,
                              'hydrogen': 0.0,
                              'oxygen': -4.40676212751828,
                              'carbon monoxide': -111262.34509634285,
                              'water': -242671.7203547155,
                              'hydrogen sulfide': -20600.0,
                              'helium': 0.0,
                              'argon': 0.0,
                              'carbon': 0.0}
                              # 'H': 218000.0,
                              # 'O': 249190.0,
                              # 'SO2': -296840.0}

        # # reactants
        # atom_list = []
        # new_dict = {}
        #
        # for key, value in comp.items():
        #     for key1, value1 in dict_components.items():
        #         if key == key1:
        #             new_dict = dict_components.get(key1, {})
        #             new_dict.update((x, y * value) for x, y in new_dict.items())
        #             atom_list.append(new_dict)
        #
        # # counter_reactants = Counter()
        # reactants_dict = Counter()
        # for d in atom_list:
        #     reactants_dict.update(d)

        # reactants
        atom_list = []

        for key, value in comp.items():
            for key1, value1 in dict_components.items():
                if key == key1:
                    atom_list.append(dict((x, y * value) for x, y in dict_components.get(key1, {}).items()))

        reactants_atom = Counter()
        for d in atom_list:
            reactants_atom.update(d)

        # products
        n_CO2 = reactants_atom["C"]
        n_SO2 = reactants_atom["S"]
        n_H2O = reactants_atom["H"] / 2
        products_dict = {'carbon dioxide': n_CO2, 'sulfur dioxide': n_SO2, 'water': n_H2O}

        # oxygen for complete combustion
        n_O = n_CO2 * 2 + n_SO2 * 2 + n_H2O * 1  # 2 is number of O atoms in CO2 AND SO2 and 1 is number of O atoms in H2O
        n_O2 = n_O / 2
        reactants_dict = deepcopy(comp)
        reactants_dict.update({'oxygen': n_O2})

        # LHV calculation
        LHV = 0
        for key, value in dict_enthalpy_mole.items():
            for key1, value1 in reactants_dict.items():
                if key == key1:
                    LHV += value * value1
            for key2, value2 in products_dict.items():
                if key == key2:
                    LHV -= value * value2

        # 298 K
        hw_liq = -285825.0
        hw_gas = -241820.0

        # 273 K
        # hw_liq = -287654.96084928664
        # hw_gas = -242628.01574091613

        HHV = LHV + (hw_gas - hw_liq) * products_dict["water"]

        if parameter == 'mass':
            # returns heating value in kJ/kg
            if hhv:
                heating_value = HHV / self.MolarMass
            else:
                heating_value = LHV / self.MolarMass
        else:
            # returns heating value in kJ/m3
            if hhv:
                heating_value = HHV * self.MolarDensity
            else:
                heating_value = LHV * self.MolarDensity

        return heating_value

    @staticmethod
    def CovertCompositionGERG(composition):
        gerg_composition = [0.0] * 22

        for gas_spice, composition in composition.items():
            gerg_composition[gerg_gas_spices.index(gas_spice)] = composition

        return gerg_composition

    def MolarMassGERG(self):
        """

        :param x:   Composition (mole fraction)
                    Do not send mole percents or mass fractions in the x() array, otherwise the output will be incorrect.
                    The sum of the compositions in the x() array must be equal to one.
                    The order of the fluids in this array is given at the top of this module.
        :return: Mm:  Molar mass (g/mol)
        """
        Mm = 0
        for i in range(1, NcGERG+1):
            Mm += self.x[i] * MMiGERG[i]
        return Mm

    def PressureGERG(self, D):
        """
        Sub PressureGERG(T, D, x, P, Z)

        Calculate pressure as a function of temperature and density.  The derivative d(P)/d(D) is also calculated
        for use in the iterative DensityGERG subroutine (and is only returned as a common variable).

        :return:        P: Pressure (kPa)
                        Z: Compressibility factor
                        dPdDsave - d(P)/d(D) [kPa/(mol/l)] (at constant temperature)
        //          - This variable is cached in the common variables for use in the iterative density solver, but not returned as an argument.
        """
        T = self.T  # Temperature (K)

        ar = self.AlpharGERG(itau=0, idelta=0, D=D)

        Z = 1 + ar[0][1]
        P = D * RGERG * T * Z
        dPdDsave = RGERG * T * (1 + 2 * ar[0][1] + ar[0][2])
        return P, Z, dPdDsave

    def DensityGERG(self, iFlag=0):
        """
        Sub DensityGERG(iFlag, T, P, x, D, ierr, herr)

        Calculate density as a function of temperature and pressure.  This is an iterative routine that calls PressureGERG
        to find the correct state point.  Generally only 6 iterations at most are required.
        If the iteration fails to converge, the ideal gas density and an error message are returned.
        No checks are made to determine the phase boundary, which would have guaranteed that the output is in the gas phase (or liquid phase when iFlag=2).
        It is up to the user to locate the phase boundary, and thus identify the phase of the T and P inputs.
        If the state point is 2-phase, the output density will represent a metastable state.

                        (An initial guess for the density can be sent in D as the negative of the guess for roots that are in the liquid phase instead of using iFlag=2)

        :param iFlag:   Set to 0 for strict pressure solver in the gas phase without checks (fastest mode, but output state may not be stable single phase)
                        Set to 1 to make checks for possible 2-phase states (result may still not be stable single phase, but many unstable states will be identified)
                        Set to 2 to search for liquid phase (and make the same checks when iFlag=1)
        :return:        D: Density (mol/l)
                        For the liquid phase, an initial value can be sent to the routine to avoid
                        a solution in the metastable or gas phases.
                        The initial value should be sent as a negative number.
        :return:        ierr: Error number (0 indicates no error)
        :return:        herr: Error message if ierr is not equal to zero
        :return:        D, ierr, herr
                        ierr: Error number (0 indicates no error)
                        herr: Error message if ierr is not equal to zero
        """

        P = self.P
        T = self.T
        x = self.x
        D = 0  # initial estimate of the density

        dPdD = 0.0
        d2PdTD = 0.0
        Cv = 0.0
        Cp = 0.0
        W = 0.0
        PP = 0.0

        ierr = 0
        herr = ""
        nFail = 0
        iFail = 0
        if self.P < epsilon:
            D = 0
            return
        tolr = 0.0000001
        Tcx, Dcx = self.PseudoCriticalPointGERG()

        if D > - epsilon:
            D = self.P / RGERG / self.T                # Ideal gas estimate for vapor phase
            if iFlag == 2:
                D = Dcx*3    # Initial estimate for liquid phase

        else:
            D = abs(D)                  # If D<0, then use as initial estimate

        plog = math.log(self.P)
        vlog = -math.log(D)
        for it in range(1, 51):
            if (vlog < -7) or (vlog > 100) or (it == 20) or (it == 30) or (it == 40) or (iFail == 1):
                # Current state is bad or iteration is taking too long.  Restart with completely different initial state
                iFail = 0
                if nFail > 2:
                    # Iteration failed (above loop did not find a solution or checks made below indicate possible 2-phase state)
                    ierr = 1
                    herr = "Calculation failed to converge in GERG method, ideal gas density returned."
                    D = P / RGERG / T
                nFail += 1
                if nFail == 1:
                    D = Dcx * 3  # If vapor phase search fails, look for root in liquid region
                elif nFail == 2:
                    D = Dcx * 2.5  # If liquid phase search fails, look for root between liquid and critical regions
                elif nFail == 3:
                    D = Dcx * 2  # If search fails, look for root in critical region

                vlog = -math.log(D)
            D = math.exp(-vlog)
            P2, Z, dPdDsave = self.PressureGERG(D)
            if (dPdDsave < epsilon) or (P2 < epsilon):
                # Current state is 2-phase, try locating a different state that is single phase
                vinc = 0.1
                if D > Dcx:
                    vinc = -0.1
                if it > 5:
                    vinc = vinc / 2
                if (it > 10) and (it < 20):
                    vinc = vinc / 5
                vlog += vinc
            else:
                # Find the next density with a first order Newton's type iterative scheme, with
                # log(P) as the known variable and log(v) as the unknown property.
                # See AGA 8 publication for further information.
                dpdlv = -D * dPdDsave     # d(p)/d[log(v)]
                vdiff = (math.log(P2) - plog) * P2 / dpdlv
                vlog += - vdiff
                if abs(vdiff) < tolr:
                    # Check to see if state is possibly 2-phase, and if so restart
                    if dPdDsave < 0:
                        iFail = 1
                    else:
                        D = math.exp(-vlog)

                        # If requested, check to see if point is possibly 2-phase
                        if iFlag > 0:
                            self.PropertiesGERG()
                            if ((PP <= 0) or (dPdD <= 0) or (d2PdTD <= 0)) or ((Cv <= 0) or (Cp <= 0) or (W <= 0)):
                                # Iteration failed (above loop did find a solution or checks made below indicate possible 2-phase state)
                                ierr = 1
                                herr = "Calculation failed to converge in GERG method, ideal gas density returned."
                                D = P / RGERG / T
                            return ierr, herr, D
                        return ierr, herr, D              # Iteration converged
        # Iteration failed (above loop did not find a solution or checks made below indicate possible 2-phase state)
        ierr = 1
        herr = "Calculation failed to converge in GERG method, ideal gas density returned."
        D = P / RGERG / T
        return ierr, herr, D

    def PropertiesGERG(self):
        """
        Sub PropertiesGERG(T, D, x, P, Z, dPdD, d2PdD2, d2PdTD, dPdT, U, H, S, Cv, Cp, W, G, JT, Kappa, A)

        Calculate thermodynamic properties as a function of temperature and density.  Calls are made to the subroutines
        ReducingParametersGERG, IdealGERG, and ResidualGERG.  If the density is not known, call subroutine DENSITY first
        with the known values of pressure and temperature.

        :param T: Temperature (K)
        :param D: Density (mol/l)
        :param x: Composition (mole fraction)
        :return:
                     P - Pressure (kPa)
                     Z - Compressibility factor
                  dPdD - First derivative of pressure with respect to density at constant temperature [kPa/(mol/l)]
                d2PdD2 - Second derivative of pressure with respect to density at constant temperature [kPa/(mol/l)^2]
                d2PdTD - Second derivative of pressure with respect to temperature and density [kPa/(mol/l)/K]
                  dPdT - First derivative of pressure with respect to temperature at constant density (kPa/K)
                     U - Internal energy (J/mol)
                     H - Enthalpy (J/mol)
                     S - Entropy [J/(mol-K)]
                    Cv - Isochoric heat capacity [J/(mol-K)]
                    Cp - Isobaric heat capacity [J/(mol-K)]
                     W - Speed of sound (m/s)
                     G - Gibbs energy (J/mol)
                    JT - Joule-Thomson coefficient (K/Pa)
                 Kappa - Isentropic Exponent
                     A - Helmholtz energy (J/mol)
        """
        # a0 = np.zeros(2+1)
        # ar = np.zeros((3+1, 3+1))

        T = self.T
        x = self.x

        # Calculate molar mass
        molar_mass = self.MolarMassGERG()

        # Calculate the density
        ierr, herr, D = self.DensityGERG()
        self.MolarDensity = D

        # Calculate the ideal gas Helmholtz energy, and its first and second derivatives with respect to temperature.
        a0 = self.Alpha0GERG()

        # Calculate the real gas Helmholtz energy, and its derivatives with respect to temperature and/or density.
        ar = self.AlpharGERG(itau=1, idelta=0, D=D)

        R = RGERG
        RT = R * T
        Z = 1 + ar[0][1]
        P = D * RT * Z
        dPdD = RT * (1 + 2 * ar[0][1] + ar[0][2])
        dPdT = D * R * (1 + ar[0][1] - ar[1][1])
        d2PdTD = R * (1 + 2 * ar[0][1] + ar[0][2] - 2 * ar[1][1] - ar[1][2])
        A = RT * (a0[0] + ar[0][0])
        G = RT * (1 + ar[0][1] + a0[0] + ar[0][0])
        U = RT * (a0[1] + ar[1][0])
        H = RT * (1 + ar[0][1] + a0[1] + ar[1][0])
        S = R * (a0[1] + ar[1][0] - a0[0] - ar[0][0])
        Cv = -R * (a0[2] + ar[2][0])
        if D > epsilon:
            Cp = Cv + T * (dPdT / D) * (dPdT / D) / dPdD
            d2PdD2 = RT * (2 * ar[0][1] + 4 * ar[0][2] + ar[0][3]) / D
            JT = (T / D * dPdT / dPdD - 1) / Cp / D  #  '=(dB/dT*T-B)/Cp for an ideal gas, but dB/dT is not known
        else:
            Cp = Cv + R
            d2PdD2 = 0
            JT = 1E+20
        W = 1000 * Cp / Cv * dPdD / molar_mass
        if W < 0:
            W = 0
        W = math.sqrt(W)
        Kappa = pow(W, 2) * molar_mass / (RT * 1000 * Z)

        self.MolarMass = molar_mass
        self.MolarDensity = D
        self.Z = Z
        self.dPdD = dPdD
        self.d2PdD2 = d2PdD2
        self.dPdT = dPdT
        self.energy = U
        self.enthalpy = H
        self.entropy = S
        self.Cv_molar = Cv
        self.Cp_molar = Cp
        self.Cv = Cv * 1000 / molar_mass
        self.Cp = Cp * 1000 / molar_mass
        self.c = W  # speed of sound [m/s]
        self.gibbs_energy = G  # Gibbs energy [J/mol]
        self.JT = JT / 1e3  # Joule-Thomson coefficient [K/Pa]
        self.isentropic_exponent = Kappa  # Isentropic exponent
        self.rho = self.MolarMass * self.P / Z /RT  # density [kg/m3]
        self.SG = self.MolarMass / air_molar_mass
        self.R_specific = R / self.MolarMass

    def Alpha0GERG(self):
        """
        Private Sub Alpha0GERG(T, D, x, a0)

        Calculate the ideal gas Helmholtz energy and its derivatives with respect to tau and delta.
        This routine is not needed when only P (or Z) is calculated.
        :param T: Temperature (K)
        :param D: Density (mol/l)
        :param x: Composition (mole fraction)
        :return: a0:        a0(0) - Ideal gas Helmholtz energy (all dimensionless [i.e., divided by RT])
                            a0(1) - tau*partial(a0)/partial(tau)
                            a0(2) - tau^2*partial^2(a0)/partial(tau)^2
        """
        T = self.T
        D = self.MolarDensity
        x = self.x
        th0T = 0.0
        LogxD = 0.0
        SumHyp0 = 0.0
        SumHyp1 = 0.0
        SumHyp2 = 0.0
        hcn = 0.0
        hsn = 0.0

        a0 = [0] * 3
        if D > epsilon:
            LogD = math.log(D)
        else:
            LogD = math.log(epsilon)
        LogT = math.log(T)
        for i in range(1, NcGERG + 1):
            if x[i] > epsilon:
                LogxD = LogD + math.log(x[i])
                SumHyp0 = 0
                SumHyp1 = 0
                SumHyp2 = 0
            for j in range(4, 8):
                if th0i[i][j] > epsilon:
                    th0T = th0i[i][j] / T
                    ep = math.exp(th0T)
                    em = 1 / ep
                    hsn = (ep - em) / 2
                    hcn = (ep + em) / 2
                    if j == 4 or j == 6:
                        LogHyp = math.log(abs(hsn))
                        SumHyp0 = SumHyp0 + n0i[i][j] * LogHyp
                        SumHyp1 = SumHyp1 + n0i[i][j] * th0T * hcn / hsn
                        SumHyp2 = SumHyp2 + n0i[i][j] * (th0T / hsn) * (th0T / hsn)
                    else:
                        LogHyp = math.log(abs(hcn))
                        SumHyp0 = SumHyp0 - n0i[i][j] * LogHyp
                        SumHyp1 = SumHyp1 - n0i[i][j] * th0T * hsn / hcn
                        SumHyp2 = SumHyp2 + n0i[i][j] * (th0T / hcn) * (th0T / hcn)

            a0[0] += +x[i] * (LogxD + n0i[i][1] + n0i[i][2] / T - n0i[i][3] * LogT + SumHyp0)
            a0[1] += +x[i] * (n0i[i][3] + n0i[i][2] / T + SumHyp1)
            a0[2] += -x[i] * (n0i[i][3] + SumHyp2)
        return a0

    # The following routines are low-level routines that should not be called outside of this code.
    def ReducingParametersGERG(self):
        """

        :param Tr:
        :param Dr:
        :return:
        // Private Sub ReducingParametersGERG(x, Tr, Dr)

        // Calculate reducing variables.  Only need to call this if the composition has changed.

        // Inputs:
        //    x() - Composition (mole fraction)

        // Outputs:
        //     Tr - Reducing temperature (K)
        //     Dr - Reducing density (mol/l)
        """
        global Drold, Trold

        # Check to see if a component fraction has changed.  If x is the same as the previous call, then exit.
        icheck = 0
        for i in range(1, NcGERG+1):
            if abs(self.x[i] - xold[i]) > 0.0000001:
                icheck = 1
            xold[i] = self.x[i]
        if icheck == 0:
            Dr = Drold
            Tr = Trold
            return Tr, Dr

        # Calculate reducing variables for T and D
        Dr = 0
        Vr = 0
        Tr = 0
        for i in range(1, NcGERG+1):
            if self.x[i] > epsilon:
                F = 1
                for j in range(i, NcGERG+1):
                    if self.x[j] > epsilon:
                        xij = F * (self.x[i] * self.x[j]) * (self.x[i] + self.x[j])
                        Vr = Vr + xij * gvij[i][j] / (bvij[i][j] * self.x[i] + self.x[j])
                        Tr = Tr + xij * gtij[i][j] / (btij[i][j] * self.x[i] + self.x[j])
                        F = 2
        if Vr > epsilon:
            Dr = 1 / Vr
        Drold = Dr
        Trold = Tr

        return Tr, Dr

    def AlpharGERG(self, itau, idelta, D):
        """
        Private Sub AlpharGERG(itau, idelta, T, D, x, ar)

        Calculate dimensionless residual Helmholtz energy and its derivatives with respect to tau and delta.

        :param itau:   Set this to 1 to calculate "ar" derivatives with respect to tau [i.e., ar(1,0), ar(1,1), and ar(2,0)], otherwise set it to 0.
        :param idelta: Currently not used, but kept as an input for future use in specifying the highest density derivative needed.
        :param T:      Temperature (K)
        :param D:      Density (mol/l)
        :param x:      Composition (mole fraction)
        :return:        ar(0,0) - Residual Helmholtz energy (dimensionless, =a/RT)
                        ar(0,1) -     delta*partial  (ar)/partial(delta)
                        ar(0,2) -   delta^2*partial^2(ar)/partial(delta)^2
                        ar(0,3) -   delta^3*partial^3(ar)/partial(delta)^3
                        ar(1,0) -       tau*partial  (ar)/partial(tau)
                        ar(1,1) - tau*delta*partial^2(ar)/partial(tau)/partial(delta)
                        ar(2,0) -     tau^2*partial^2(ar)/partial(tau)^2
        """

        T = self.T
        x = self.x

        global Told, Trold, Trold2, Drold

        global Tr, Dr
        delp = [0] * (7+1)
        Expd = [0] * (7+1)
        ar = [[0] * 4 for _ in range(4)]

        for i in range(4):
            for j in range(4):
                ar[i][j] = 0

        # Set up del, tau, log(tau), and the first 7 calculations for del^i
        Tr, Dr = self.ReducingParametersGERG()
        delta = D / Dr
        tau = Tr / T
        lntau = math.log(tau)
        delp[1] = delta
        Expd[1] = math.exp(-delp[1])
        for i in range(2, 8):
            delp[i] = delp[i - 1] * delta
            Expd[i] = math.exp(-delp[i])

        # If temperature has changed, calculate temperature dependent parts
        if (abs(T - Told) > 0.0000001) or (abs(Tr - Trold2) > 0.0000001):
            self.tTermsGERG(lntau, x)
        Told = T
        Trold2 = Tr

        # Calculate pure fluid contributions
        for i in range(1, NcGERG+1):
            if x[i] > epsilon:
                for k in range(1, int(kpol[i]+1)):
                    ndt = x[i] * delp[int(doik[i][k])] * taup[i][k]
                    ndtd = ndt * doik[i][k]
                    ar[0][1] += ndtd
                    ar[0][2] += ndtd * (doik[i][k] - 1)
                    if itau > 0:
                        ndtt = ndt * toik[i][k]
                        ar[0][0] += ndt
                        ar[1][0] += ndtt
                        ar[2][0] += ndtt * (toik[i][k] - 1)
                        ar[1][1] += ndtt * doik[i][k]
                        ar[1][2] += ndtt * doik[i][k] * (doik[i][k] - 1)
                        ar[0][3] += ndtd * (doik[i][k] - 1) * (doik[i][k] - 2)
                for k in range(int(kpol[i]+1), int(kpol[i] + kexp[i]+1)):
                    ndt = x[i] * delp[int(doik[i][k])] * taup[i][k]*Expd[int(coik[i][k])]
                    ex = coik[i][k] * delp[int(coik[i][k])]
                    ex2 = doik[i][k] - ex
                    ex3 = ex2 * (ex2 - 1)
                    ar[0][1] += ndt * ex2
                    ar[0][2] += ndt * (ex3 - coik[i][k] * ex)
                    if itau > 0:
                        ndtt = ndt * toik[i][k]
                        ar[0][0] += ndt
                        ar[1][0] += ndtt
                        ar[2][0] += ndtt * (toik[i][k] - 1)
                        ar[1][1] += ndtt * ex2
                        ar[1][2] += ndtt * (ex3 - coik[i][k] * ex)
                        ar[0][3] += ndt * (ex3 * (ex2 - 2) - ex * (3 * ex2 - 3 + coik[i][k]) * coik[i][k])

        # Calculate mixture contributions
        for i in range(1, NcGERG):  # for (int i = 1; i <= NcGERG - 1; ++i)
            if x[i] > epsilon:
                for j in range(i+1, NcGERG+1):  # for (int j = i + 1; j <= NcGERG; ++j)
                    if x[j] > epsilon:
                        mn = int(mNumb[i][j])
                        if mn >= 0:
                            xijf = x[i] * x[j] * fij[i][j]
                            for k in range(1, int(kpolij[mn] + 1)):  # for (int k = 1; k <= kpolij[mn]; ++k)
                                ndt = xijf * delp[int(dijk[mn][k])] * taupijk[mn][k]
                                ndtd = ndt * dijk[mn][k]
                                ar[0][1] += ndtd
                                ar[0][2] += ndtd * (dijk[mn][k] - 1)
                                if itau > 0:
                                    ndtt = ndt * tijk[mn][k]
                                    ar[0][0] += ndt
                                    ar[1][0] += ndtt
                                    ar[2][0] += ndtt * (tijk[mn][k] - 1)
                                    ar[1][1] += ndtt * dijk[mn][k]
                                    ar[1][2] += ndtt * dijk[mn][k] * (dijk[mn][k] - 1)
                                    ar[0][3] += ndtd * (dijk[mn][k] - 1) * (dijk[mn][k] - 2)

                            for k in range(int(1+kpolij[mn]), int(kpolij[mn]+kexpij[mn]+1)):  # for (int k = 1 + kpolij[mn]; k <= kpolij[mn] + kexpij[mn]; ++k)
                                cij0 = cijk[mn][k] * delp[2]
                                eij0 = eijk[mn][k] * delta
                                ndt = xijf * nijk[mn][k] * delp[int(dijk[mn][k])] * math.exp(cij0 + eij0 + gijk[mn][k] + tijk[mn][k] * lntau)
                                ex = dijk[mn][k] + 2 * cij0 + eij0
                                ex2 = (ex * ex - dijk[mn][k] + 2 * cij0)
                                ar[0][1] += ndt * ex
                                ar[0][2] += ndt * ex2
                                if itau > 0:
                                    ndtt = ndt * tijk[mn][k]
                                    ar[0][0] += ndt
                                    ar[1][0] += ndtt
                                    ar[2][0] += ndtt * (tijk[mn][k] - 1)
                                    ar[1][1] += ndtt * ex
                                    ar[1][2] += ndtt * ex2
                                    ar[0][3] += ndt * (ex * (ex2 - 2 * (dijk[mn][k] - 2 * cij0)) + 2 * dijk[mn][k])
        return ar


    def tTermsGERG(self, lntau, x):
        """
        // Private Sub tTermsGERG(lntau, x)

        // Calculate temperature dependent parts of the GERG-2008 equation of state
        :param lntau:
        :param x:
        :return:
        """
        taup0 = [0] * (12+1)

        i = 5  # Use propane to get exponents for short form of EOS
        for k in range(1, int(kpol[i] + kexp[i] + 1)):  # for (int k = 1; k <= kpol[i] + kexp[i]; ++k)
            taup0[k] = math.exp(toik[i][k] * lntau)
        for i in range(1, NcGERG+1):  # for (int i = 1; i <= NcGERG; ++i)
            if x[i] > epsilon:
                if (i > 4) and (i != 15) and (i != 18) and (i != 20):
                    for k in range(1, int(kpol[i]+kexp[i]+1)):  # for (int k = 1; k <= kpol[i] + kexp[i]; ++k)
                        taup[i][k] = noik[i][k] * taup0[k]
                else:
                    for k in range(1, int(kpol[i]+kexp[i]+1)):  # for (int k = 1; k <= kpol[i] + kexp[i]; ++k)
                        taup[i][k] = noik[i][k] * math.exp(toik[i][k] * lntau)

        for i in range(1, NcGERG):  # for (int i = 1; i <= NcGERG - 1; ++i)
            if x[i] > epsilon:
                for j in range(i+1, NcGERG+1):  # for (int j = i + 1; j <= NcGERG; ++j)
                    if x[j] > epsilon:
                        mn = int(mNumb[i][j])
                        if mn >= 0:
                            for k in range(1, int(kpolij[mn]+1)):  # for (int k = 1; k <= kpolij[mn]; ++k)
                                taupijk[mn][k] = nijk[mn][k] * math.exp(tijk[mn][k] * lntau)

    def PseudoCriticalPointGERG(self):
        """
        // PseudoCriticalPointGERG(x, Tcx, Dcx)

        // Calculate a pseudo critical point as the mole fraction average of the critical temperatures and critical volumes
        :param x:
        :return:
        """
        x = self.x

        Vcx = 0
        Tcx = 0
        Vcx = 0
        Dcx = 0
        for i in range(1, NcGERG+1):  # for (int i = 1; i <= NcGERG; ++i)
            Tcx = Tcx + x[i] * Tc[i]
            Vcx = Vcx + x[i] / Dc[i]
        if Vcx > epsilon:
            Dcx = 1 / Vcx

        return Tcx, Dcx
