#   #!/usr/bin/env python
#   -*- coding: utf-8 -*-
#   ******************************************************************************
#     Copyright (c) 2022.
#     Developed by Yifei Lu
#     Last change on 7/29/22, 9:06 AM
#     Last change by yifei
#    *****************************************************************************


class GasStorage:
    def __init__(self, capacity, initial_contents=None):
        self.capacity = capacity
        self.contents = {} if initial_contents is None else initial_contents
        self.composition = self.calculate_stored_gas_composition()

    def store_gas(self, gas_composition):
        for gas_type, amount in gas_composition.items():
            print(gas_type, amount)
            if gas_type in self.contents.keys():
                self.contents[gas_type] += amount
            else:
                self.contents[gas_type] = amount

        if sum(self.contents.values()) > self.capacity:
            print("Warning: Storage capacity exceeded.")

        self.composition = self.calculate_stored_gas_composition()

    def withdraw_gas(self, amount):
        for k, v in self.composition.items():
            self.contents[k] -= amount * v
            if self.contents[k] < 0:
                self.contents[k] = 0

    def calculate_stored_gas_composition(self):
        s = sum(self.contents.values())

        composition = {}

        for k, v in self.contents.items():
            pct = v / s
            composition[k] = pct
        return composition

    def get_fill_level(self):
        return sum(self.contents.values()) / self.capacity

    def get_total_contents(self):
        return self.contents

    def get_gas_amount(self, gas_type):
        return self.contents.get(gas_type, 0)


if __name__ == "__main__":
    # Example usage
    initial_contents = {"methane": 100, "propane": 50}
    storage = GasStorage(capacity=300, initial_contents=initial_contents)

    print("Initial contents:", storage.get_total_contents())

    storage.store_gas({"butane": 75})
    storage.store_gas({"methane": 50})
    storage.withdraw_gas(30)

    print("Updated contents:", storage.get_total_contents())

    print("Amount of methane:", storage.get_gas_amount("methane"))

    print(f"Fill level of the storage: {storage.get_fill_level()}")

    print(f"Stored gas composition: {storage.calculate_stored_gas_composition()}")
