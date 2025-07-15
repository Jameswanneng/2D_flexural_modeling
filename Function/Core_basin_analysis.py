#!/usr/bin/env python3
# =======================================
# Well Data Processing Module - Handles well location and stratigraphic data

# Key Functionality:
#- Parses well location files with distance measurements
#- Reads age-dependent stratigraphic data from Excel files
#- Maintains consistent data structures for basin modeling

# (c) China University of Geosciences (Beijing)
#            ALL RIGHTS RESERVED
#        Authors: Neng Wan, Xueyan Li


import matplotlib
import re
import xlrd
import numpy as np

# Configure matplotlib for publication-quality output
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'Arial'


def get_well_locations(well_list_file):
    """
        Reads well names and their corresponding distances from a text file.

        Input:
            well_list_file: Path to the well list file. Each line contains a well name and distance.

        Returns:
            dict: A dictionary mapping well names to their corresponding distances (float).
        """
    well_list_file = open(well_list_file, 'r', encoding='utf-8')
    well_list = []
    well_distance_dict = {}
    well_location_list = []
    for line in well_list_file:
        if line.startswith('#') or line.startswith('>'):
            continue
        line = re.split('\s+', line)
        well_name = line[0]
        distance = float(line[1])
        well_distance_dict[well_name] = distance
        well_location_list.append(distance)
        well_list.append(well_name)
    well_list_file.close()
    return well_distance_dict


def read_age_dependant_data(basin_age_sequence, base_age, excel_file_path, column_index=1):
    """
    Reads age-dependent data from an Excel file based on a given age sequence.

    Input:
        basin_age_sequence (list): List of ages representing different stages of basin evolution.
        base_age (float or int): The reference age used to determine the starting point in the sequence.
        excel_file_path (str): Path to the Excel file containing the data in Sheet1.
        column_index (int): Starting column index for reading data (default is 1, i.e., second column).

    Returns:
        dict: A dictionary where keys are ages and values are NumPy arrays of corresponding column data.
    """
    workbook = xlrd.open_workbook(excel_file_path, 'r')
    data_sheet = workbook.sheet_by_name('Sheet1')
    age_dependant_data_dict = {}
    base_age_index = basin_age_sequence.index(base_age)
    for age_index_0 in range(base_age_index, len(basin_age_sequence)):
        age_0 = basin_age_sequence[age_index_0]
        age_dependant_data_dict[age_0] = np.array(data_sheet.col_values(column_index, start_rowx=1, end_rowx=None))
        column_index += 1
    return age_dependant_data_dict
