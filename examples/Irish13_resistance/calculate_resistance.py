#   #!/usr/bin/env python
#   -*- coding: utf-8 -*-
#   ******************************************************************************
#     Copyright (c) 2022.
#     Developed by Yifei Lu
#     Last change on 5/27/22, 1:44 PM
#     Last change by yifei
#    *****************************************************************************
import pandas as pd

import GasNetSim as gns
from pathlib import Path

network = gns.create_network_from_csv(Path('../Irish13/'))
network_resistance = pd.DataFrame(columns=['resistance_index', 'inlet_index', 'outlet_index', 'resistance'])

i = 0
for pipe in network.pipelines.values():
    resistance = pipe.calculate_fictitious_resistance()
    new_row = pd.DataFrame([[i+1, pipe.inlet_index, pipe.outlet_index, resistance]],
                           columns=['resistance_index', 'inlet_index', 'outlet_index', 'resistance'])
    network_resistance = pd.concat([network_resistance, new_row], ignore_index=True)
    i += 1

network_resistance.to_csv("Irish13_resistance.csv", sep=';')
