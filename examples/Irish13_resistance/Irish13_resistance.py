#   #!/usr/bin/env python
#   -*- coding: utf-8 -*-
#   ******************************************************************************
#     Copyright (c) 2022.
#     Developed by Yifei Lu
#     Last change on 6/13/22, 4:47 PM
#     Last change by yifei
#    *****************************************************************************
import GasNetSim as gns
from pathlib import Path

network = gns.create_network_from_csv(Path('.'))
network.simulation()
