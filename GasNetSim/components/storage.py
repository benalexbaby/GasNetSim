#   #!/usr/bin/env python
#   -*- coding: utf-8 -*-
#   ******************************************************************************
#     Copyright (c) 2022.
#     Developed by Yifei Lu
#     Last change on 7/29/22, 9:06 AM
#     Last change by yifei
#    *****************************************************************************
from .node import Node


class Storage(Node):
    def __init__(self, node_index):
        super().__init__(node_index)
