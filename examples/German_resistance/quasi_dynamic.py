#   #!/usr/bin/env python
#   -*- coding: utf-8 -*-
#   ******************************************************************************
#     Copyright (c) 2022.
#     Developed by Yifei Lu
#     Last change on 8/2/22, 3:32 PM
#     Last change by yifei
#    *****************************************************************************
import GasNetSim as gns
from pathlib import Path

network = gns.create_network_from_csv(Path('.'))
gns.run_time_series(network, 'profiles/test_profiles.csv')
