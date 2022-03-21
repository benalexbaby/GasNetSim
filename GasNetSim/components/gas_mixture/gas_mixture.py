#  #!/usr/bin/env python
#  -*- coding: utf-8 -*-
#  ******************************************************************************
#    Copyright (c) 2021.
#    Developed by Yifei Lu
#    Last change on 12/21/21, 4:22 PM
#    Last change by yifei
#   *****************************************************************************
from collections import OrderedDict

from .thermo.thermo import Mixture


class GasMixture:
    """
    Class for gas mixture properties
    """
    def __init__(self, pressure, temperature, composition: OrderedDict, method='PREOS'):
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
            self.thermo_Mixture = Mixture(P=pressure, T=temperature, zs=composition)

    @property
    def compressibility(self):
        if self.method == "PREOS":
            return self.thermo_Mixture.Z