#  #!/usr/bin/env python
#  -*- coding: utf-8 -*-
#  ******************************************************************************
#    Copyright (c) 2021.
#    Developed by Yifei Lu
#    Last change on 12/21/21, 4:58 PM
#    Last change by yifei
#   *****************************************************************************
from collections import OrderedDict


# NATURAL_GAS = OrderedDict([('methane', 0.96522),
#                            ('nitrogen', 0.00259),
#                            ('carbon dioxide', 0.00596),
#                            ('ethane', 0.01819),
#                            ('propane', 0.0046),
#                            ('isobutane', 0.00098),
#                            ('butane', 0.00101),
#                            ('2-methylbutane', 0.00047),
#                            ('pentane', 0.00032),
#                            ('hexane', 0.00066)])


NATURAL_GAS_gri30 = OrderedDict([('methane', 0.9477),
                                 ('ethane', 0.042),
                                 ('propane', 0.002),
                                 ('nitrogen', 0.005),
                                 ('carbon dioxide', 0.003),
                                 ('oxygen', 0.0001),
                                 ('hydrogen', 0.0002)])

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

HYDROGEN = OrderedDict([('hydrogen', 1)])
