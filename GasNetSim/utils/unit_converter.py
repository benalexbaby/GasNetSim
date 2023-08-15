#   #!/usr/bin/env python
#   -*- coding: utf-8 -*-
#   ******************************************************************************
#     Copyright (c) 2023.
#     Developed by Yifei Lu
#     Last change on 3/9/23, 4:46 PM
#     Last change by yifei
#    *****************************************************************************

def unit_converter(value: float, unit_from: str, unit_to: str):
    from pint import UnitRegistry
    ureg = UnitRegistry()
    Q_ = ureg.Quantity
    converted_value = Q_(f"{value} * {unit_from}").to(f"{unit_to}")
    return converted_value.magnitude


if __name__ == "__main__":
    print(unit_converter(value=255659,
                         unit_from="GWh/year",
                         unit_to="GWh/d"))
