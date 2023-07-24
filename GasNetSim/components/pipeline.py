#   #!/usr/bin/env python
#   -*- coding: utf-8 -*-
#   ******************************************************************************
#     Copyright (c) 2022.
#     Developed by Yifei Lu
#     Last change on 6/19/22, 11:11 AM
#     Last change by yifei
#    *****************************************************************************
import logging
import math
from scipy.constants import bar, atm, zero_Celsius

from .node import *
from .utils.pipeline_function.flow_rate import *
from .utils.pipeline_function.friction_factor import *
from .utils.pipeline_function.outlet_temperature import *
from .utils.gas_mixture.gas_mixture import *


STANDARD_TEMPERATURE = zero_Celsius + 15  # 15 degree Celsius in Kelvin
STANDARD_PRESSURE = 1 * atm  # 1 atm in pa


class Pipeline:
    """
    Class for gas transmission pipelines
    """

    def __init__(self, inlet: Node, outlet: Node, diameter, length, efficiency=0.85, roughness=0.000015,
                 ambient_temp=15+zero_Celsius, ambient_pressure=1*atm, heat_transfer_coefficient=3.69, valve=0,
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
        self.resistance = self.calculate_fictitious_resistance()
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
            if t1 == ambient_temp or t2 == ambient_temp:
                return self.inlet.temperature
            else:
                return ambient_temp + (t1 - t2) / math.log((t1 - ambient_temp) / (t2 - ambient_temp))
        except (ZeroDivisionError, ValueError):
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
        specific_gravity = self.gas_mixture.specific_gravity
        h1 = self.inlet.altitude
        h2 = self.outlet.altitude
        avg_pressure = self.calc_average_pressure()
        avg_temperature = self.calc_average_temperature()
        z = self.gas_mixture.compressibility

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
        pb = STANDARD_PRESSURE
        tb = STANDARD_TEMPERATURE
        if flow_rate is not None:
            cross_section = math.pi * (self.diameter/2)**2
            return flow_rate * pb/tb * self.calc_average_temperature() / self.calc_average_pressure()/cross_section
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
                                  rho=self.gas_mixture.density, viscosity=self.gas_mixture.viscosity)
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
        implemented_methods = ['weymouth', 'chen', 'nikuradse', 'colebrook-white', 'hagen-poiseuille']

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
        elif method == "hagen-poiseuille":
            hagen_poiseuille(N_re=self.calculate_reynolds_number())

    def calculate_fictitious_resistance(self):
        tb = STANDARD_TEMPERATURE
        pb = STANDARD_PRESSURE
        d = self.diameter
        length = self.length
        e = self.efficiency
        return pb * length**0.5 / (FLOW_EQUATION_CONSTANT * tb * d**2.5 * e)

    def calc_physical_char_gas_pipe(self):
        """
        Calculate physical characteristics of the gas pipe which is a combined term of all variables not impacted by
        the gas transmission state variables
        :return: Gas pipeline physical characteristics
        """

        tb = STANDARD_TEMPERATURE  # Temperature base, 15 Celsius
        pb = STANDARD_PRESSURE  # Pressure base, 1 atm
        d = self.diameter
        length = self.length
        # avg_temperature = self.calc_average_temperature()
        # z = self.gas_mixture.compressibility
        e = self.efficiency
        # f = self.calculate_pipe_friction_factor()
        # if f is None:
        #     f = 0.01
        # specific_gravity = self.gas_mixture.specific_gravity

        # if specific_gravity < 0:
        #     specific_gravity = 0.5
        #     logging.debug(self.gas_mixture.zs)
        #     logging.warning("Gas mixture specific gravity is smaller than 0, set it as default value 0.5.")

        return (FLOW_EQUATION_CONSTANT * tb / pb) * (d ** 2.5) * ((1 / length)**0.5) * e

    def calculate_coefficient_for_iteration(self):
        avg_temperature = self.calc_average_temperature()
        z = self.gas_mixture.compressibility
        f = self.calculate_pipe_friction_factor()
        if f is None:
            f = 0.01
        specific_gravity = self.gas_mixture.specific_gravity

        if specific_gravity < 0:
            specific_gravity = 0.5
            logging.debug(self.gas_mixture.zs)
            logging.warning("Gas mixture specific gravity is smaller than 0, set it as default value 0.5.")

        pipeline_physical_characteristic = self.calc_physical_char_gas_pipe()

        return pipeline_physical_characteristic / (specific_gravity * avg_temperature * z * f)**0.5

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
            print(f"Pipeline {self.inlet_index} has same pressure on both ends!")
            # raise ValueError('Got condition case 0.')
            return 0

    def calc_flow_rate(self):
        """
        Calculate the volumetric flow rate through the pipe
        :return: Volumetric flow rate [sm3/s]
        """
        flow_direction = self.determine_flow_direction()
        p1 = self.inlet.pressure
        p2 = self.outlet.pressure

        slope_correction = self.calc_pipe_slope_correction()
        tmp = self.calculate_coefficient_for_iteration()

        return flow_direction * abs(p1 ** 2 - p2 ** 2 - slope_correction) ** (1 / 2) * tmp

    def flow_rate_first_order_derivative(self, is_inlet=True):
        p1 = self.inlet.pressure
        p2 = self.outlet.pressure
        slope_corr = self.calc_pipe_slope_correction()
        pipeline_coefficient = self.calculate_coefficient_for_iteration()
        tmp = (abs(p1 ** 2 - p2 ** 2 - slope_corr)) ** (-0.5)

        if is_inlet:
            return pipeline_coefficient * p1 * tmp
        else:
            return pipeline_coefficient * p2 * tmp

    def calc_gas_mass_flow(self):
        """
        Calculate gas mass flow rate through the pipe
        :return: Mass flow rate [kg/s]
        """
        q = self.calc_flow_rate()
        gas_rho = GasMixture(composition=self.get_mole_fraction(),
                             pressure=STANDARD_PRESSURE,
                             temperature=STANDARD_TEMPERATURE).density
        return q * gas_rho

    def calc_pipe_outlet_temp(self):
        """
        Calculate pipe outlet temperature based on the physical law of flow temperature loss
        :return: Pipe outlet temperature
        """
        qm = abs(self.calc_gas_mass_flow())
        friction = self.calculate_pipe_friction_factor()
        if qm is not None and friction is not None and self.gas_mixture.heat_capacity_constant_pressure is not None:
            beta = calculate_beta_coefficient(ul=3.69,
                                              qm=qm,
                                              cp=self.gas_mixture.heat_capacity_constant_pressure,
                                              D=self.diameter)
            gamma = calculate_gamma_coefficient(mu_jt=self.gas_mixture.joule_thomson_coefficient,
                                                Z=self.gas_mixture.compressibility,
                                                R_specific=self.gas_mixture.R_specific,
                                                f=friction,
                                                qm=qm,
                                                p_avg=self.calc_average_pressure(),
                                                D=self.diameter)
            return calculate_pipeline_outlet_temperature(beta=beta,
                                                         gamma=gamma,
                                                         Ts=self.ambient_temp,
                                                         L=self.length,
                                                         T1=self.inlet.temperature)
        else:
            return self.ambient_temp

    def get_mole_fraction(self):
        """
        Get mole fraction of the gas composition inside pipeline
        :return: Gas mole fraction
        """
        # mole_fraction = dict()
        # for i in range(len(self.gas_mixture.IDs)):
        #     gas = self.gas_mixture.IDs[i]
        #     try:
        #         mole_fraction[gas] = self.gas_mixture.zs[i]
        #     except TypeError:
        #         print(mole_fraction)
        # return mole_fraction
        return self.gas_mixture.composition


class Resistance:
    def __init__(self, inlet: Node, outlet: Node, resistance=1e6):
        self.inlet = inlet
        self.outlet = outlet
        self.inlet_index = inlet.index
        self.outlet_index = outlet.index
        self.resistance = resistance
        self.flow_rate = None
        self.gas_mixture = self.inlet.gas_mixture

    def update_gas_mixture(self):
        self.gas_mixture = self.inlet.gas_mixture

    def calc_pipe_slope_correction(self):
        return 0

    def calculate_coefficient_for_iteration(self):
        avg_temperature = 288.15  # TODO temperature calculation for resistance
        z = self.gas_mixture.compressibility
        f = 0.01
        specific_gravity = self.gas_mixture.specific_gravity

        if specific_gravity < 0:
            specific_gravity = 0.5
            logging.debug(self.gas_mixture.zs)
            logging.warning("Gas mixture specific gravity is smaller than 0, set it as default value 0.5.")

        return 1 / (specific_gravity * avg_temperature * z * f)**0.5 / self.resistance

    def determine_flow_direction(self):
        """
        Determine the flow direction inside a pipeline
        :return: -1 or 1, respectively from inlet to outlet or contrariwise
        """
        p1 = self.inlet.pressure
        p2 = self.outlet.pressure
        slope_correction = 0  # TODO slope correction equals 0?

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

        slope_correction = 0
        tmp = self.calculate_coefficient_for_iteration()

        return flow_direction * abs(p1 ** 2 - p2 ** 2 - slope_correction) ** (1 / 2) * tmp

    def flow_rate_first_order_derivative(self, is_inlet=True):
        p1 = self.inlet.pressure
        p2 = self.outlet.pressure
        slope_corr = self.calc_pipe_slope_correction()
        pipeline_coefficient = self.calculate_coefficient_for_iteration()
        tmp = (abs(p1 ** 2 - p2 ** 2 - slope_corr)) ** (-0.5)

        if is_inlet:
            return pipeline_coefficient * p1 * tmp
        else:
            return pipeline_coefficient * p2 * tmp


    def calc_gas_mass_flow(self):
        """
        Calculate gas mass flow rate through the pipe
        :return: Mass flow rate [kg/s]
        """
        q = self.calc_flow_rate()
        gas_rho = GasMixture(composition=self.gas_mixture.composition,
                             pressure=STANDARD_PRESSURE,
                             temperature=STANDARD_TEMPERATURE).density
        return q * gas_rho

    def calc_pipe_outlet_temp(self):
        """

        :return:
        """
        # TODO: non-static method
        return 288.15


class LinearResistance:
    def __init__(self, inlet: Node, outlet: Node, resistance=1e6):
        self.inlet = inlet
        self.outlet = outlet
        self.inlet_index = inlet.index
        self.outlet_index = outlet.index
        self.resistance = resistance
        self.flow_rate = None
        self.gas_mixture = self.inlet.gas_mixture

    def update_gas_mixture(self):
        self.gas_mixture = self.inlet.gas_mixture

    def calc_pipe_slope_correction(self):
        return 0

    def calculate_coefficient_for_iteration(self):
        avg_temperature = 288.15  # TODO temperature calculation for resistance
        z = self.gas_mixture.compressibility
        f = 0.01
        specific_gravity = self.gas_mixture.specific_gravity

        if specific_gravity < 0:
            specific_gravity = 0.5
            logging.debug(self.gas_mixture.zs)
            logging.warning("Gas mixture specific gravity is smaller than 0, set it as default value 0.5.")

        return 1 / (specific_gravity * avg_temperature * z * f)**0.5 / self.resistance

    def determine_flow_direction(self):
        """
        Determine the flow direction inside a pipeline
        :return: -1 or 1, respectively from inlet to outlet or contrariwise
        """
        p1 = self.inlet.pressure
        p2 = self.outlet.pressure

        try:
            p1 - p2
        except ValueError or TypeError:
            print(f'p1: {p1}, p2: {p2}')
        if p1 > p2:
            return 1
        elif p1 < p2:
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

        tmp = self.calculate_coefficient_for_iteration()

        return flow_direction * abs(p1 - p2) * tmp

    def flow_rate_first_order_derivative(self, is_inlet=True):
        p1 = self.inlet.pressure
        p2 = self.outlet.pressure
        pipeline_coefficient = self.calculate_coefficient_for_iteration()

        return pipeline_coefficient


    def calc_gas_mass_flow(self):
        """
        Calculate gas mass flow rate through the pipe
        :return: Mass flow rate [kg/s]
        """
        q = self.calc_flow_rate()
        gas_rho = GasMixture(composition=self.gas_mixture.composition,
                             pressure=STANDARD_PRESSURE,
                             temperature=STANDARD_TEMPERATURE).density
        return q * gas_rho

    def calc_pipe_outlet_temp(self):
        """

        :return:
        """
        # TODO: non-static method
        return 288.15


class ShortPipe:
    def __init__(self, inlet: Node, outlet: Node):
        self.inlet = inlet
        self.outlet = outlet
        self.inlet_index = inlet.index
        self.outlet_index = outlet.index
        self.flow_rate = - self.inlet.volumetric_flow  # Short pipes are used to connect to supply nodes
        self.gas_mixture = self.inlet.gas_mixture

    def update_gas_mixture(self):
        self.gas_mixture = self.inlet.gas_mixture

    def calc_pipe_slope_correction(self):
        return 0

    def calculate_coefficient_for_iteration(self):
        return 0

    def determine_flow_direction(self):
        """
        Determine the flow direction inside a pipeline
        :return: -1 or 1, respectively from inlet to outlet or contrariwise
        """
        p1 = self.inlet.pressure
        p2 = self.outlet.pressure
        slope_correction = 0  # TODO slope correction equals 0?

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

        return - self.inlet.volumetric_flow

    def calc_gas_mass_flow(self):
        """
        Calculate gas mass flow rate through the pipe
        :return: Mass flow rate [kg/s]
        """
        q = self.calc_flow_rate()
        gas_rho = GasMixture(composition=self.gas_mixture.composition,
                             pressure=STANDARD_PRESSURE,
                             temperature=STANDARD_TEMPERATURE).density
        return q * gas_rho

    def calc_pipe_outlet_temp(self):
        """

        :return:
        """
        # TODO: non-static method
        return 288.15
