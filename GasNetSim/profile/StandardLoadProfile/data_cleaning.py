#  #!/usr/bin/env python
#  -*- coding: utf-8 -*-
#  ******************************************************************************
#    Copyright (c) 2022.
#    Developed by Yifei Lu
#    Last change on 3/7/22, 11:48 AM
#    Last change by yifei
#   *****************************************************************************
import pandas as pd
import datetime
import os
from pathlib import Path
import glob

# Define constants to be used
N_H_CLASS = 11  # Number of hausehold classes
N_G_CLASS = 4  # Number of businuss building classes
N_TEMP_CASES = 10  # Number of temperature ranges
N_WIND_CASES = 5


def pct_dist_gas_consumption(filename):
    # Read csv-file
    df = pd.read_csv(filename, sep=';', encoding='ISO-8859-1')

    # Convert all percentage strings into float numbers
    for col in df.columns[2:]:
        try:
            df[col] = df[col].apply(lambda x: x.replace(',', '.'))
            df[col] = df[col].apply(lambda x: float(x[:-1]) / 100)
        except AttributeError:
            # print(df[col])
            print("Data in " + filename.split("\\")[-1] + " are not string type.")
            print(col)

    df = df[df.columns[2:]]

    try:
        df.columns = [int(x.split('.')[0]) for x in df.columns]
    except ValueError:
        print(filename.split("\\")[-1])
    cols = list(range(24))
    df = df[cols]

    df.columns = pd.Series(df.columns).apply(lambda x: datetime.datetime.strptime(str(x) + ":30", '%H:%M'))

    # Set a time range for interpolation
    times = pd.date_range(start='1/1/1900 00:30:00', end='1/2/1900 00:30:00', freq='15min')

    extended_df = pd.DataFrame(columns=list(times))

    for col in extended_df.columns:
        if col in df.columns:
            extended_df[col] = df[col].copy()

    extended_df[extended_df.columns[-1]] = extended_df[extended_df.columns[0]]

    # Change all columns into numeric types, so that the interpolation is able to be done
    for col in extended_df:
        extended_df[col] = pd.to_numeric(extended_df[col], errors='coerce')

    # Interpolate data with an order of 3
    extended_df = extended_df.interpolate(method='polynomial', order=3, axis=1)
    extended_df = extended_df.drop(extended_df.columns[-1], axis=1)

    extended_df.columns = [x - datetime.timedelta(days=1) if x.day == 2 else x for x in extended_df.columns]

    # reorder the clock hours of a gas day time into natural day
    times = pd.date_range(start='1/1/1900 00:00:00', end='1/1/1900 23:59:00', freq='15min')
    extended_df = extended_df[times]

    extended_df = extended_df.div(4.0)

    summary_dfs = dict()
    index_list = list()

    if "familienhaeuser" in filename:

        for n in range(N_H_CLASS):
            index_list.append([x + N_TEMP_CASES * n for x in list(range(N_TEMP_CASES))])

        class_names = ["Class_" + str(x + 1) for x in range(N_H_CLASS)]
        for n in range(N_H_CLASS):
            summary_dfs[class_names[n]] = extended_df.iloc[index_list[n], :].reset_index(drop=True).copy()

    elif "Gewerbe" in filename:
        for n in range(N_G_CLASS):
            index_list.append([x + N_TEMP_CASES * n for x in list(range(N_TEMP_CASES))])

        class_names = ["Class_" + str(x + 1) for x in range(N_G_CLASS)]
        for n in range(N_G_CLASS):
            summary_dfs[class_names[n]] = extended_df.iloc[index_list[n], :].reset_index(drop=True).copy()
    else:
        print("Check filename: " + filename)

    return summary_dfs


def hh_slp():
    dir_path = os.path.abspath('')
    data_dir = Path(dir_path + "/data")
    filenames = glob.glob(str(data_dir) + "/*.csv")

    summary = dict()

    for f in filenames:
        if "familienhaeuser" in f:
            class_to_slp = pct_dist_gas_consumption(f)
            summary[f.split('\\')[-1].split('.')[0]] = class_to_slp

    return summary


def gewerbe_slp():
    dir_path = os.path.abspath('')
    data_dir = Path(dir_path + "/data")
    filenames = glob.glob(str(data_dir) + "/*.csv")

    summary = dict()

    for f in filenames:
        if "Gewerbe" in f:
            class_to_slp = pct_dist_gas_consumption(f)
            summary[f.split('\\')[-1].split('.')[0]] = class_to_slp

    return summary


def temp_classification(amp_temp=8):
    class_type = 0

    if amp_temp <= -15:
        class_type = 0
    elif -15 < amp_temp <= -10:
        class_type = 1
    elif -10 < amp_temp <= -5:
        class_type = 2
    elif -5 < amp_temp <= 0:
        class_type = 3
    elif 0 < amp_temp <= 5:
        class_type = 4
    elif 5 < amp_temp <= 10:
        class_type = 5
    elif 10 < amp_temp <= 15:
        class_type = 6
    elif 15 < amp_temp <= 20:
        class_type = 7
    elif 20 < amp_temp <= 25:
        class_type = 8
    elif amp_temp > 25:
        class_type = 9
    else:
        print("Please check the temperature input!")

    return class_type


def sigmoid_coef(filename):
    # read csv-file
    df = pd.read_csv(Path('data/EFH_sigmoid.csv'), sep=';', encoding='ISO-8859-1')

    # reset column names
    df.columns = df.iloc[0, :]
    df = df.drop(0).reset_index(drop=True)

    # Convert elements in DataFrame to float type
    for col in df.columns[2:]:
        try:
            df[col] = df[col].apply(lambda x: x.replace(',', '.'))
            df[col] = df[col].apply(lambda x: float(x))
        except AttributeError:
            # print(df[col])
            print("Data in " + filename.split("\\")[-1] + " are not string type.")
            print(col)

    df = df[df.columns[2:]]

    index_list = list()

    for n in range(N_H_CLASS):
        index_list.append([x + N_WIND_CASES * n for x in list(range(N_WIND_CASES))])

    summary_dfs = dict()

    class_names = ["Class_" + str(x + 1) for x in range(N_H_CLASS)]

    for n in range(N_H_CLASS):
        summary_dfs[class_names[n]] = df.iloc[index_list[n], :].reset_index(drop=True).copy()

    return summary_dfs


if __name__ == "__main__":
    hh_slp_dict = hh_slp()
    gewerbe_slp_dict = gewerbe_slp()

    ambient_temp = float(input("Set the ambient temperature: "))
    print(temp_classification(ambient_temp))
