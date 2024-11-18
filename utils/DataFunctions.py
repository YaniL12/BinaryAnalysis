from astropy.io import fits
import numpy as np
import pandas as pd
import os
import re

def FitsToDF(fn):
    """
    Converts a FITS file to a pandas DataFrame.

    Parameters:
    fn (str): The path to the FITS file.

    Returns:
    pandas.DataFrame: The DataFrame containing the data from the FITS file.
    """
    hdul = fits.open(fn)
    data = hdul[1].data
    return pd.DataFrame(data)

def FitsToDFWithVariableLengthCols(fn):
    """
    Converts a FITS file with columns of varying lengths to a pandas DataFrame.

    Parameters:
    fn (str): The path to the FITS file.

    Returns:
    tuple: A tuple containing the DataFrame and a dictionary of columns with varying lengths.
    """
    with fits.open(fn) as hdul:
        data = hdul[1].data
        data_dict = {}
        variable_length_cols = {}
        
        for name in data.names:
            col_data = data[name]
            if isinstance(col_data[0], (np.ndarray, list)):
                variable_length_cols[name] = col_data
            else:
                data_dict[name] = col_data
        
        df = pd.DataFrame(data_dict)
        
        return df, variable_length_cols


def parse_readme(dat_file, readme_file):
    """
    Parses the readme file to extract column specifications and names.

    Parameters:
    dat_file (str): The path to the data file.
    readme_file (str): The name of the readme file.

    Returns:
    tuple: A tuple containing the column specifications and column names.
    """
    readme_file = os.path.dirname(dat_file) + '/' + readme_file

    with open(readme_file, 'r') as file:
        lines = file.readlines()

    table4_found = False
    byte_description = []
    data_line_start = np.nan
    table_name = os.path.basename(dat_file)
    first_line_text = 'Byte-by-byte Description of file: ' + table_name

    table_cols = []

    # Find the start of the byte-by-byte description for table4.dat
    for i, line in enumerate(lines):
        if first_line_text in line:
            table4_found = True

            # Get the table column names. These should be 2 lines after the first_line_text
            table_cols = lines[i + 2].split()
        
        if table4_found:       
            if line.__contains__('1-'):
                data_line_start = i
                break

    for i, line in enumerate(lines[data_line_start:]):
        if line.startswith('----'):
            break
        else:
            byte_description.append(line.strip())

    colspecs = []
    column_label = []
    start = 0


    label_index = table_cols.index('Label')

    for line in byte_description:
        parts = line.split()
        if '-' in parts[0]:
            # This is the byte range. E.g. 1 - 8
            start = int(parts[0].split('-')[0])

            # These are columns where the byte range has gone from space separated to hyphen separated with no spaces. E.g. 1 - 8 to 1-8:
            if len(parts[0].split('-')[1]) > 0:
                end = int(parts[0].split('-')[1])
                column_label.append(parts[label_index ])
            else:
                end = int(parts[1])
                column_label.append(parts[label_index + 1])

            colspecs.append((start - 1, end))
            # column_label.append(' '.join(parts[label_index + 1:]))

            # print(column_label[-1])

        # What if the first element is a number? E.g. the first byte index is just = 5 not 0-5
        elif parts[0].isdigit() and '-' in parts:
            end = int(parts[0])
            colspecs.append((start - 1, end))
            column_label.append(parts[label_index + 1])

        else:
            continue
    
    return colspecs, column_label

def read_dat_file(dat_file, readme_file="ReadMe"):
    """
    Reads a data file and returns a DataFrame based on the column specifications in the readme file.

    Parameters:
    dat_file (str): The path to the data file.
    readme_file (str): The name of the readme file. Default is "ReadMe".

    Returns:
    pandas.DataFrame: The DataFrame containing the data from the data file.
    """
    try:
        colspecs, column_names = parse_readme(dat_file, readme_file)
        df = pd.read_fwf(dat_file, colspecs=colspecs, names=column_names)
        return df
    except pd.errors.EmptyDataError:
        print("Error: The file is empty or no columns to parse.")
    except Exception as e:
        print(f"An error occurred: {e}")



# Define a custom function to handle the data parsing
def custom_split(line):
    # First, split the line by commas
    parts = re.split(r',\s*', line)
    
    # For the second element, further split by spaces
    if len(parts) > 1:
        second_element_split = parts[1].split()
        # Replace the second element with the two parts split by space
        parts = parts[:1] + second_element_split + parts[2:]
    
    return parts

def read_binary_result_file(fn, coltype=0):
    """
    Reads a custom result file and returns a DataFrame.

    Parameters:
    fn (str): The path to the custom result file.

    Returns:
    pandas.DataFrame: The DataFrame containing the data from the custom result file.
    """

    # Read the data from the text file
    with open(fn, 'r') as file:
        data = [custom_split(line.strip()) for line in file]


# 'f_contr': 0.5, 'mass_1': 0.888822, 'rv_1': 116.17699, 'vmic_1': 1.5, 'vsini_1': 4.0, 'mass_2': 0.888822, 'rv_2': 38.0, 'vmic_2': 1.5, 'vsini_2': 4.0, 'teff_1': 5.3533408203125, 'teff_2': 5.3533408203125, 'logg_1': 4.327199, 'logg_2': 4.327199, 'logl_1': 0, 'logl_2': 0}

    if coltype == 0:
        cols = [
            'sobject_id',
            'residual',
            'rchi2',
            'f_contr',
            'mass_1',
            'age_1',
            # 'metallicity_1',
            'rv_1',
            'fe_h_1',
            'vmic_1',
            'vsini_1',
            'mass_2',
            'age_2',
            # 'metallicity_2',
            'rv_2',
            'fe_h_2',
            'vmic_2',
            'vsini_2',
            'teff_1',
            'teff_2',
            'logg_1',
            'logg_2',
            'logl_1',
            'logl_2'
        ]
    elif coltype == 1:
        cols = [
            'sobject_id',
            'residual',
            'rchi2',
            'f_contr',
            'mass_1',
            'rv_1',
            'vmic_1',
            'vsini_1',
            'mass_2',
            'rv_2',
            'vmic_2',
            'vsini_2',
            'teff_1',
            'teff_2',
            'logg_1',
            'logg_2',
            'logl_1',
            'logl_2',
            'age_1',
            'age_2',
            'FeH_1',
            'FeH_2'
        ]
    else:
        cols = ['sobject_id', 'residual', 'rchi2', 'f_contr', 'mass_1', 'age_1', 'metallicity_1', 'rv_1', 'fe_h_1', 'vmic_1', 'vsini_1', 'teff_1', 'logg_1', 'logl_1', 'mass_2', 'age_2', 'metallicity_2', 'rv_2', 'fe_h_2', 'vmic_2', 'vsini_2', 'teff_2', 'logg_2', 'logl_2']

    # Convert the data to a pandas DataFrame
    data = pd.DataFrame(data, columns=cols)
    # Remove rows with None or NaN in 'age_1' column
    data = data.dropna(subset=['age_1'])
    # Convert all columns except the first one to float
    data.iloc[:, 0] = data.iloc[:, 0].astype(int)
    data.iloc[:, 1:] = data.iloc[:, 1:].astype(float)
    data['delta_rv_GALAH'] = abs(data['rv_2'] - data['rv_1'])

    return data
