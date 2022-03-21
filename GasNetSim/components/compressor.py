#  #!/usr/bin/env python
#  -*- coding: utf-8 -*-
#  ******************************************************************************
#    Copyright (c) 2022.
#    Developed by Yifei Lu
#    Last change on 1/17/22, 11:21 AM
#    Last change by yifei
#   *****************************************************************************
from .node import *


class Compressor:
    """
    Class to formulate compressor stations.
    """
    def __init__(self, inlet: Node, outlet: Node,
                 drive='electric',
                 compression_ratio=1.1,
                 thermodynamic_process='isentropic'):
        self.cp = inlet.gas_mixture.cp
        self.cv = inlet.gas_mixture.cv
        self.T1 = inlet.temperature
        self.r_comp = compression_ratio
        if thermodynamic_process == 'isentropic':
            self.n = self.cp / self.cv
        elif thermodynamic_process == 'isothermal':
            self.n = 1
        else:
            raise ValueError('Only isentropic or isothermal process is currently supported.')

    def power_consumption(self):
        return self.cp*self.T1*(self.r_comp**((self.n-1)/self.n) - 1)*self.qm