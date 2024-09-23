from astropy.io import fits
import numpy as np
import pandas as pd
import os

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

    # Find the start of the byte-by-byte description for table4.dat
    for i, line in enumerate(lines):
        if first_line_text in line:
            table4_found = True
        
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
    column_names = []
    start = 0

    for line in byte_description:
        parts = line.split()
        if '-' in parts[0]:
            start = int(parts[0].split('-')[0])
            end = int(parts[1])

            colspecs.append((start - 1, end))
            column_names.append(parts[4])

        elif parts[0].isdigit():
            end = int(parts[0])
            colspecs.append((start - 1, end))
            column_names.append(parts[3])

        else:
            continue
    
    return colspecs, column_names

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