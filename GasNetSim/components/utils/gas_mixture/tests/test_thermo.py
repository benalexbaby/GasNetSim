#   #!/usr/bin/env python
#   -*- coding: utf-8 -*-
#   ******************************************************************************
#     Copyright (c) 2022.
#     Developed by Yifei Lu
#     Last change on 12/20/22, 10:35 AM
#     Last change by yifei
#    *****************************************************************************
import thermo

print(thermo.__version__)

ONE_ATM = 101325
ITERATIONS = 100

results = list()
for i in range(ITERATIONS + 1):
    composition = {"methane": 1 - i/ITERATIONS, "hydrogen": i / ITERATIONS}
    gas_mix = thermo.Mixture(T=300, P=50*ONE_ATM, zs=composition)
    gas_mix.atomss
    gas_mix.set_eos()
    results.append(gas_mix.Z)
    # plt.plot(results)