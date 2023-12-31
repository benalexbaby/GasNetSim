#  #!/usr/bin/env python
#  -*- coding: utf-8 -*-
#  ******************************************************************************
#    Copyright (c) 2021.
#    Developed by Yifei Lu
#    Last change on 12/21/21, 4:22 PM
#    Last change by yifei
#   *****************************************************************************
from collections import OrderedDict
import logging

from .thermo.thermo import Mixture
# from thermo import Mixture


class GasMixture:
    """
    Class for gas mixture properties
    """
    def __init__(self, pressure, temperature, composition, method='PREOS'):
        """

        :param pressure:
        :param temperature:
        :param composition:
        :param method:
        """
        self.pressure = pressure
        self.temperature = temperature
        self.composition = composition
        self.method = method
        if method == "GERG-2008":
            pass
        elif method == "PREOS":
            self.thermo_mixture = Mixture(P=pressure, T=temperature, zs=composition)

    @property
    def compressibility(self):
        if self.method == "PREOS":
            if self.thermo_mixture.Z is not None:
                z = self.thermo_mixture.Z
            else:
                logging.warning("Compressibility is not available, using the Z for gas!")
                z = self.thermo_mixture.Zg
            return self.thermo_mixture.Z

    @property
    def specific_gravity(self):
        if self.method == "PREOS":
            if self.thermo_mixture.SG is not None:
                specific_gravity = self.thermo_mixture.SG
            else:
                logging.warning("Specific gravity is not available, using the SG for gas!")
                specific_gravity = self.thermo_mixture.SGg
            return self.thermo_mixture.SG

    @property
    def molar_mass(self):
        if self.method == "PREOS":
            return self.thermo_mixture.MW

    @property
    def density(self):
        if self.method == "PREOS":
            return self.thermo_mixture.rho

    @property
    def joule_thomson_coefficient(self):
        if self.method == "PREOS":
            return self.thermo_mixture.JT

    @property
    def viscosity(self):
        if self.method == "PREOS":
            return self.thermo_mixture.mu

    @property
    def heat_capacity_constant_pressure(self):
        if self.method == "PREOS":
            return self.thermo_mixture.Cp

    @property
    def R_specific(self):
        if self.method == "PREOS":
            return self.thermo_mixture.R_specific


