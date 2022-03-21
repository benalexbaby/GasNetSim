# import sys
# import os
# sys.path.append(os.path.dirname(__file__))

import GasNetSim as gns
from pathlib import Path

network = gns.create_network_from_csv(Path('.'))
network.simulation()
