#  #!/usr/bin/env python
#  -*- coding: utf-8 -*-
#  ******************************************************************************
#    Copyright (c) 2022.
#    Developed by Yifei Lu
#    Last change on 3/21/22, 3:22 PM
#    Last change by yifei
#   *****************************************************************************
import math
import logging


def calculate_height_difference_correction(h1, h2, z, p, t, sg):
    if h1 == h2:
        return 0

    try:
        return 0.06835 * sg * (h2 - h1) * p ** 2 / (z * t)
    except:
        print("Error calculating the slope correction!")
