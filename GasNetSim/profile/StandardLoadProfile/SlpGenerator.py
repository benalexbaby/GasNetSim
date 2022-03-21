#  #!/usr/bin/env python
#  -*- coding: utf-8 -*-
#  ******************************************************************************
#    Copyright (c) 2022.
#    Developed by Yifei Lu
#    Last change on 3/7/22, 11:48 AM
#    Last change by yifei
#   *****************************************************************************
import pandas as pd
import numpy as np
from data_cleaning import *
from matplotlib import pyplot as plt
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()


def calc_temp(Td, Td_1, Td_2, Td_3):
    """
    Moving average temperature
    :param Td:
    :param Td_1:
    :param Td_2:
    :param Td_3:
    :return:
    """
    return (Td + 0.5 * Td_1 + 0.25 * Td_2 + 0.125 * Td_3) / (1 + 0.5 + 0.25 + 0.125)


def h_value(A, B, C, D, theta, theta_0=40.0, m_h=0, m_w=0, b_h=0, b_w=0):
    """
    Corrected value for daily demand
    :param A:
    :param B:
    :param C:
    :param D:
    :param theta:
    :param theta_0:
    :param m_h:
    :param m_w:
    :param b_h:
    :param b_w:
    :return:
    """
    f_sigmoid = A/(1+(B/(theta-theta_0))**C) + D
    f_linear = max((m_h*theta + b_h), (m_w*theta + b_w))
    return f_sigmoid + f_linear


def daily_demand(kw, h, f=1):
    """

    :param kw: Kundenswert
    :param h: h_value
    :param f: Wochentagsfaktoren (around 1)
    :return:
    """
    return kw*h*f


def customer_value(Q_array, theta_array):
    h_sum = 0
    for theta in theta_array:
        h_sum += h_value(3.1850191, -37.4124155, 6.1723179, 0.0761096, theta, 40.0, 0, 0, 0, 0)
    return Q_array.sum()/h_sum


def average_daily_ambient_temp(file, region='DEA2', frequency="hour"):
    """

    :param file: file path to temperature file
    :param region: NUTS2 identifier of a city or region. default DEA2 (Aachen region)
    :param frequency: time interval between data entires
    :return:
    """
    df = pd.read_csv(file)
    filtered_cols = ['utc_timestamp']

    for col in df.columns:
        if region in col:
            filtered_cols.append(col)

    regional_weather_data = df[filtered_cols]

    if frequency == "day":

        daily_temp = list()

        dates = pd.date_range(start='2014-12-28', end='2015-12-31', freq='D')

        dates = dates.to_list()

        for i in range(len(dates)):
            daily_temp.append(regional_weather_data[region + '_temperature'][i*24:i*24 + 23].sum()/24)

        df_daily_temp = pd.DataFrame({'dates': dates, 'temp': daily_temp})

        return df_daily_temp
    elif frequency == "hour":
        df_hourly_temp = pd.DataFrame({'time': regional_weather_data['utc_timestamp'],
                                       'temp': regional_weather_data['DEA2_temperature']})
        return df_hourly_temp


if __name__ == '__main__':
    T_allokation = calc_temp(-2.0, 0.5, 3.4, 3.6)
    print("AllokationsTemperatur ist: " + str(T_allokation))
    h = h_value(3.1850191, -37.4124155, 6.1723179, 0.0761096, T_allokation, 40.0, 0, 0, 0, 0)
    print("h-Wert ist: " + str(h))
    Q = daily_demand(50, h, 1)
    print("Tagesmenge ist: " + str(Q))

    temp_dea2 = average_daily_ambient_temp("weather_data_filtered.csv")
    sigmoid_coeffients = sigmoid_coef("data/EFH_sigmoid.csv")

    nrw_sigmoid_coef_efh = sigmoid_coeffients['Class_3'].iloc[1, :]

    gas_consumption = pd.DataFrame()

    kundenwerte = [200, 100, 300, 300, 400, 150]
    gas_consumption['time'] = range(1, 366)

    for j in range(len(kundenwerte)):
        tmp = list()
        for i in range(4, 369):
            T_ma = calc_temp(temp_dea2['temp'][i], temp_dea2['temp'][i-1],
                                     temp_dea2['temp'][i-2], temp_dea2['temp'][i-3])
            h = h_value(nrw_sigmoid_coef_efh[0], nrw_sigmoid_coef_efh[1], nrw_sigmoid_coef_efh[2], nrw_sigmoid_coef_efh[3],
                        T_ma, 40.0, 0, 0, 0, 0)
            Q = daily_demand(kundenwerte[j], h, 1)
            tmp.append(Q)
        gas_consumption[f'{j+2}'] = tmp

    plt.figure()
    plt.plot(gas_consumption)
    plt.show()

    gas_consumption.to_csv('gas_profile.csv')


    # plt.plot(temp_dea2['dates'], temp_dea2['temp'])

    # gewerbe_slp_dict = gewerbe_slp()
    # for key, value in gewerbe_slp_dict.items():
    #     temp_class = temp_classification(T_allokation)
    #     result_slp = Q * value['Class_4'].iloc[temp_class, :]
    #     plt.figure()
    #     plt.plot(result_slp)
    #     plt.show()

