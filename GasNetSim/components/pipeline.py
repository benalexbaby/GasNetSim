#  #!/usr/bin/env python
#  -*- coding: utf-8 -*-
#  ******************************************************************************
#    Copyright (c) 2021.
#    Developed by Yifei Lu
#    Last change on 12/21/21, 4:17 PM
#    Last change by yifei
#   *****************************************************************************
import logging
import math

from .node import Node
from .pipeline_function.friction_factor import *
from .pipeline_function.outlet_temperature import *
from .gas_mixture.thermo.thermo import Mixture


class Pipeline:
    """
    Class for gas transmission pipelines
    """

    def __init__(self, inlet: Node, outlet: Node, diameter, length, efficiency=0.85, roughness=0.000015,
                 ambient_temp=288.15, ambient_pressure=101325, heat_transfer_coefficient=3.69, valve=0,
                 friction_factor_method='chen'):
        """

        :param inlet: Gas pipe inlet node
        :param outlet: Gas pipe outlet node
        :param diameter: Pipe diameter [m]
        :param length: Pipe length [m]
        :param efficiency: Pipe transmission efficiency, default 0.85 (normal conditioned)
        :param ambient_temp: Pipe surrounding temperature [K]
        :param ambient_pressure: Pipe surrounding temperature [Pa]
        """

        self.inlet = inlet
        self.outlet = outlet
        self.inlet_index = inlet.index
        self.outlet_index = outlet.index
        self.diameter = diameter
        self.length = length
        self.efficiency = efficiency
        self.ambient_temp = ambient_temp
        self.ambient_pressure = ambient_pressure
        self.flow_rate = None
        self.mass_flow_rate = None
        self.roughness = roughness
        self.valve = valve
        self.gas_mixture = self.inlet.gas_mixture
        self.friction_factor_method = friction_factor_method

    def update_gas_mixture(self):
        self.gas_mixture = self.inlet.gas_mixture

    def calc_average_temperature(self):
        """
        Calculate average gas temperature inside pipe
        :return: Average gas temperature [K]
        """

        try:
            ambient_temp = self.ambient_temp
            t1 = self.inlet.temperature
            t2 = self.outlet.temperature
            return ambient_temp + (t1 - t2) / math.log((t1 - ambient_temp) / (t2 - ambient_temp))
        except ZeroDivisionError:
            return self.inlet.temperature
        except ValueError:
            return self.inlet.temperature

    def calc_average_pressure(self):
        """
        Calculate average gas pressure inside pipe
        :return: Average gas pressure
        """
        try:
            p1 = self.inlet.pressure
        except TypeError:
            return None
        try:
            p2 = self.outlet.pressure
        except TypeError:
            return None

        return 2./3. * ((p1 + p2) - (p1 * p2) / (p1 + p2))


    def calc_pipe_slope_correction(self):
        """
        Calculate the slope correction factor, which caused by the inclination of the pipe and adds up to the effect
        caused by the pressure difference
        :return: Slope correction factor
        """

        if self.gas_mixture.SG is not None:
            specific_gravity = self.gas_mixture.SG
        else:
            logging.warning("Specific gravity is not available, using the SG for gas!")
            specific_gravity = self.gas_mixture.SGg
        h1 = self.inlet.altitude
        h2 = self.outlet.altitude
        avg_pressure = self.calc_average_pressure()
        avg_temperature = self.calc_average_temperature()
        if self.gas_mixture.Z is not None:
            z = self.gas_mixture.Z
        else:
            z = self.gas_mixture.Zg

        try:
            return 0.06835 * specific_gravity * (h2 - h1) * avg_pressure ** 2 / (z * avg_temperature)
        except:
            print("Error calculating the slope correction!")

    def calc_flow_velocity(self):
        """
        Calculate flow velocity in m/s
        :return:
        """
        flow_rate = self.flow_rate
        if flow_rate is not None:
            cross_section = math.pi * (self.diameter/2)**2
            return flow_rate * 101325 / 288.15 * self.calc_average_temperature() / self.calc_average_pressure() / cross_section
        else:
            return None

    def calculate_reynolds_number(self):
        """
        Calculate the Reynolds number
        :return: Reynolds number
        """
        if self.flow_rate is not None:
            flow_velocity = self.calc_flow_velocity()
            return reynold_number(diameter=self.diameter, velocity=flow_velocity,
                                  rho=self.gas_mixture.rho, viscosity=self.gas_mixture.mu)
        else:
            # if the flow rate cannot be calculated yet, set the Reynolds number to be 1e7
            return 1e7

    def calculate_pipe_friction_factor(self):
        """
        Calculate pipe friction factor, inside fully turbulent flow field the friction factor can be simplified only
        related to pipe diameter. Some works choose a fix number as 0.01.
        :return: Pipe friction factor
        """
        # Friction factor models implemented in this tool, if not available, add new methods/models in the
        # friction_factor.py file
        implemented_methods = ['weymouth', 'chen', 'nikuradse', 'colebrook-white']

        method = self.friction_factor_method

        if method in implemented_methods:
            pass
        else:
            raise ValueError(f'Friction calculation method {method} is not defined! Choose on from '
                             f'{implemented_methods} or implement you own in friction_factor.py')

        if method == 'weymouth':
            return 0.0093902 / (self.diameter ** (1 / 3))
        elif method == 'chen':
            return chen(epsilon=self.roughness, d=self.diameter, N_re=self.calculate_reynolds_number())
        elif method == 'nikuradse':
            return nikuradse(d=self.diameter, epsilon=self.roughness)
        elif method == 'colebrook-white':
            return colebrook_white(epsilon=self.roughness, d=self.diameter, N_re=self.calculate_reynolds_number())

    def calc_physical_char_gas_pipe(self):
        """
        Calculate physical characteristics of the gas pipe which is a combined term of all variables not impacted by
        the gas transmission state variables
        :return: Gas pipeline physical characteristics
        """

        tb = 15 + 273.15  # Temperature base, 15 Celsius
        pb = 101325  # Pressure base, 1 bar
        d = self.diameter
        length = self.length
        avg_temperature = self.calc_average_temperature()
        if self.gas_mixture.Z is not None:
            z = self.gas_mixture.Z
        else:
            z = self.gas_mixture.Zg
        e = self.efficiency
        f = self.calculate_pipe_friction_factor()
        if f is None:
            f = 0.01
        if self.gas_mixture.SG is not None:
            specific_gravity = self.gas_mixture.SG
        else:
            specific_gravity = self.gas_mixture.SGg

        if specific_gravity < 0:
            specific_gravity = 0.5
            print(self.gas_mixture.zs)
            print("Gas mixture specific gravity is smaller than 0, set it as default value 0.5.")

        return (13.29 * tb / pb) * (d ** 2.5) * \
               ((1 / (length * specific_gravity * avg_temperature * z * f)) ** 0.5) * e

    def determine_flow_direction(self):
        """
        Determine the flow direction inside a pipeline
        :return: -1 or 1, respectively from inlet to outlet or contrariwise
        """
        p1 = self.inlet.pressure
        p2 = self.outlet.pressure
        slope_correction = self.calc_pipe_slope_correction()
        try:
            p1 ** 2 - p2 ** 2 - slope_correction
        except ValueError or TypeError:
            print(f'p1: {p1}, p2: {p2}')
        if p1 ** 2 - p2 ** 2 - slope_correction > 0:
            return 1
        elif p1 ** 2 - p2 ** 2 - slope_correction < 0:
            return -1
        else:
            raise ValueError('Got condition case 0.')

    def calc_flow_rate(self):
        """
        Calculate the volumetric flow rate through the pipe
        :return: Volumetric flow rate [sm3/s]
        """
        flow_direction = self.determine_flow_direction()
        p1 = self.inlet.pressure
        p2 = self.outlet.pressure
        slope_correction = self.calc_pipe_slope_correction()
        pipe_physical_char = self.calc_physical_char_gas_pipe()

        return flow_direction * abs(p1 ** 2 - p2 ** 2 - slope_correction) ** (1 / 2) * pipe_physical_char

    def calc_gas_mass_flow(self):
        """
        Calculate gas mass flow rate through the pipe
        :return: Mass flow rate [kg/s]
        """
        q = self.calc_flow_rate()
        gas_rho = Mixture(zs=self.get_mole_fraction(), P=101325, T=288.15).rho
        return q * gas_rho

    def calc_pipe_outlet_temp(self):
        """
        Calculate pipe outlet temperature based on the physical law of flow temperature loss
        :return: Pipe outlet temperature
        """
        qm = self.calc_gas_mass_flow()
        friction = self.calculate_pipe_friction_factor()
        if qm is not None and friction is not None and self.gas_mixture.Cp is not None:
            beta = calc_beta(ul=3.69, qm=qm, cp=self.gas_mixture.Cp, d=self.diameter)
            gamma = calc_gamma(mu_jt=self.gas_mixture.JT, z=self.gas_mixture.Z, R=self.gas_mixture.R_specific,
                               f=friction, qm=qm, p_a=self.calc_average_pressure(), D=self.diameter)
            return outlet_temp(beta=beta, gamma=gamma, Ts=self.ambient_temp, L=self.length, T1=self.inlet.temperature)
        else:
            return self.ambient_temp


    def get_mole_fraction(self):
        """
        Get mole fraction of the gas composition inside pipeline
        :return: Gas mole fraction
        """
        mole_fraction = dict()
        for i in range(len(self.gas_mixture.components)):
            gas = self.gas_mixture.components[i]
            try:
                mole_fraction[gas] = self.gas_mixture.zs[i]
            except TypeError:
                print(mole_fraction)
        return mole_fraction
