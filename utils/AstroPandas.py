# Converts FITS files to Pandas DataFrames


from astropy.io import fits
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

def FitsToDF(fn):
    hdul = fits.open(fn)
    data = hdul[1].data
    # Correct for 'endian' data. See https://stackoverflow.com/questions/30283836/creating-pandas-dataframe-from-numpy-array-leads-to-strange-errors

    return pd.DataFrame((np.array(data).byteswap().newbyteorder()))


    try:
        df = pd.DataFrame(data)
        return pd.DataFrame(df)
    except Exception as e:
        print("Big endian - little endian error. Correcting...")
        return pd.DataFrame((np.array(data).byteswap().newbyteorder()))

def FitsToDFWithVariableLengthCols(fn):
    with fits.open(fn) as hdul:
        data = hdul[1].data
        data_dict = {}
        variable_length_cols = {}
        
        for name in data.names:
            col_data = data[name]
            if isinstance(col_data[0], (np.ndarray, list)):
                # Handle columns with arrays of varying lengths
                variable_length_cols[name] = col_data
                # print(f"Handling column {name} with varying lengths separately")
            else:
                data_dict[name] = col_data
        
        df = pd.DataFrame(data_dict)
        
        return df, variable_length_cols
    

def PlotSpectra(data_or_wavelengths, flux=None, wr=None):
    x = np.empty
    y= np.empty
    
    if isinstance(data_or_wavelengths, pd.DataFrame):
        x, y = data_or_wavelengths['wave'].values, data_or_wavelengths['sob'].values  # Convert to numpy arrays
    elif data_or_wavelengths is not None and flux is not None:
        x, y = np.array(data_or_wavelengths), np.array(flux)  # Convert to numpy arrays
    else:
        raise ValueError("Invalid arguments. Provide either a DataFrame or two lists/arrays.")
    
    
    plt.figure(figsize=(25,5))

    if wr is not None:
        mask = (x >= wr[0]) & (x <= wr[1])
        plt.plot(x[mask], y[mask])
    else:
        plt.plot(x, y)

    plt.xlabel('Wavelength (A)')
    plt.ylabel('Flux')
    plt.show()


# To plot data from multiple CCDs, ignoring the gaps between them.
def split_data_by_gaps(wavelengths, flux, gap_threshold=5):
    # Find the indices where the difference between consecutive wavelengths exceeds the gap_threshold
    gaps = np.where(np.diff(wavelengths) > gap_threshold)[0]
    
    # Split the data at the gaps
    split_indices = np.split(np.arange(len(wavelengths)), gaps + 1)
    
    return [(wavelengths[indices], flux[indices]) for indices in split_indices]
