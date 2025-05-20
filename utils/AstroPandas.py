# Converts FITS files to Pandas DataFrames


from astropy.io import fits
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

import sys 
import os

working_directory = '/avatar/yanilach/PhD-Home/binaries_galah-main/spectrum_analysis/BinaryAnalysis'
sys.path.append(os.path.join(working_directory, 'utils'))
import DataFunctions as df
from astropy.table import Table

GALAH_DR4 = None
isochrone_data = None

def FitsToDF(fn, data_index=None):
    hdul = fits.open(fn)
    if data_index is not None:
        data = hdul[data_index].data
    else:
        try:
            data = hdul[0].data
        except IndexError:
            data = hdul[1].data

    # Correct for 'endian' data. See https://stackoverflow.com/questions/30283836/creating-pandas-dataframe-from-numpy-array-leads-to-strange-errors
    # If using numpy > 2.0.0, use this line instead:
    if np.__version__ > '2.0.0':
        return pd.DataFrame(np.array(data).byteswap().view(np.array(data).dtype.newbyteorder()))
    else:
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
        
        # Correct for 'endian' data. See https://stackoverflow.com/questions/30283836/creating-pandas-dataframe-from-numpy-array-leads-to-strange-errors
        # If using numpy > 2.0.0, use this line instead:
        if np.__version__ > '2.0.0':
            df = pd.DataFrame(np.array(data).byteswap().view(np.array(data).dtype.newbyteorder()))
        else:
            df = pd.DataFrame((np.array(data).byteswap().newbyteorder()))

        return df, variable_length_cols

def LoadSpectraFitsFile(fn):

    hdul = fits.open(fn)

    try :
        data = hdul[1].data  # Assuming the spectrum is in the first extension
        header = hdul[1].header
    except IndexError:
        data = hdul[0].data  # Assuming the spectrum is in the primary HDU at index 0
        header = hdul[0].header

    start_wavelength = header['CRVAL1']  # Starting wavelength
    step = header['CDELT1']  # Wavelength step per pixel
    num_pixels = data.shape[0]  # Number of data points

    wavelength = start_wavelength + np.arange(num_pixels) * step

    if np.__version__ > '2.0.0':
        wavelength = np.array(wavelength).byteswap().view(np.array(wavelength).dtype.newbyteorder())
        data =   np.array(data).byteswap().view(np.array(data).dtype.newbyteorder())
    else:
        wavelength = (np.array(wavelength).byteswap().newbyteorder())
        data = (np.array(data).byteswap().newbyteorder())

    data =  [wavelength, data]
    data = np.array(data).T

    return pd.DataFrame(data, columns=['wave', 'sob'], index=None)  
    

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



def plot_hr_GALAH(id=None, binary_data=None, lines=False, error_map=False, components=False):

    global GALAH_DR4
    global isochrone_data

    if type(GALAH_DR4) == type(None):
        print("Loading GALAH DR4 data into cache...")
        GALAH_DR4_dir = '/avatar/buder/GALAH_DR4/'
        GALAH_DR4 = df.FitsToDF(GALAH_DR4_dir + "catalogs/galah_dr4_allspec_240207.fits")

    if type(isochrone_data) == type(None):
        print("Loading isochrone data into cache...")
        isochrone_data = Table.read('/avatar/yanilach/PhD-Home/binaries_galah-main/spectrum_analysis/BinaryAnalysis/' + 'assets/' + 'parsec_isochrones_reduced.fits')

    GALAH_data = GALAH_DR4

    plt.figure(figsize=(15, 7))
    cleaned_data = GALAH_data[['teff', 'logg']].replace([np.inf, -np.inf], np.nan).dropna()

    # Plot the 2D histogram
    plt.hexbin(cleaned_data['teff'], cleaned_data['logg'], gridsize=100, cmap='gray_r', alpha=0.5)

    age = 9.3
    isochrone_data_selection = isochrone_data[(isochrone_data['logAge'] <= age + 0.01) & (isochrone_data['logAge'] >= age - 0.01) & (isochrone_data['m_h'] == 0.0)]



    ms_binaries = [(3000, 4.4), (3500, 4.4), (4000, 4.4), (4500, 4.4), (5000, 4.3), (5500, 4.2), (5750, 4)]
    min_teff = np.min([x for x, y in ms_binaries])
    max_teff = np.max([x for x, y in ms_binaries])

    # Polyfit these coordinates
    x, y = zip(*ms_binaries)
    coefficients = np.polyfit(x, y, 2)
    polynomial = np.poly1d(coefficients)
    # Plot the polynomial
    plt.plot(np.arange(3500, max_teff), polynomial(np.arange(3500, max_teff)), color='red', label='Binary Main Sequence', ls='--')


    plt.plot(10 ** isochrone_data_selection['logT'], isochrone_data_selection['logg'], alpha=0.5, color='black', label='MS', ls='dotted')

    # Select the sample stars in the galah data and show in blue
    if id is not None:
        sample = GALAH_data[GALAH_data['sobject_id'] == id]
        plt.scatter(sample['teff'], sample['logg'], alpha=0.5, color='red', s=25, label='Unresolved')

        if binary_data is not None:
            if binary_data[binary_data['sobject_id'] == id].shape[0] > 0:
                sample = binary_data[binary_data['sobject_id'] == id]
                plt.scatter(sample['teff_1'] * 1e3, sample['logg_1'], alpha=0.5, color='blue', s=25, label='Primary')
                plt.scatter(sample['teff_2'] * 1e3, sample['logg_2'], alpha=0.5, color='green', s=25, label='Secondary')
            else:
                print("No binary data found for sobject_id", id)

    elif id is None and binary_data is not None and error_map is False:
        # Plot everything with the binary stars
        sample = GALAH_data[GALAH_data['sobject_id'].isin(binary_data['sobject_id'])]

        if components is False:
            plt.scatter(sample['teff'], sample['logg'], alpha=0.5, color='grey', s=3, label='Unresolved')
            plt.scatter(binary_data['teff_1'] * 1e3, binary_data['logg_1'], alpha=0.5, color='lightskyblue', s=5, label='Primary')
            plt.scatter(binary_data['teff_2'] * 1e3, binary_data['logg_2'], alpha=0.5, color='coral', s=5, label='Secondary')
        elif components == 'primary':
            plt.scatter(binary_data['teff_1'] * 1e3, binary_data['logg_1'], alpha=0.5, color='coral', s=5, label='Primary')
        elif components == 'secondary':
            plt.scatter(binary_data['teff_2'] * 1e3, binary_data['logg_2'], alpha=0.5, color='mediumseagreen', s=5, label='Secondary')

        if lines:
            for i in range(len(sample)):
                plt.arrow(
                    sample['teff'].iloc[i], sample['logg'].iloc[i], 
                    binary_data['teff_1'].iloc[i] * 1e3 - sample['teff'].iloc[i],  # Fix here
                    binary_data['logg_1'].iloc[i] - sample['logg'].iloc[i], 
                    color='blue', alpha=0.1, ls='dashed'
                )
                plt.arrow(
                    sample['teff'].iloc[i], sample['logg'].iloc[i], 
                    binary_data['teff_2'].iloc[i] * 1e3 - sample['teff'].iloc[i],  # Fix here
                    binary_data['logg_2'].iloc[i] - sample['logg'].iloc[i], 
                    color='green', alpha=0.1, ls='dashed'
                )

    elif error_map is not None:
        # Get a combine dataframe containing all the data from the GALAH DR4 and the binary data
        combined_data = GALAH_data.merge(binary_data, how='inner', on='sobject_id')
        combined_data = combined_data.dropna(subset=['rchi2'])
        # print(error_map)
        if error_map is not True:
            # If error_map is a string, use it as the column name
            # Check if the column exists in the combined data
            if type(error_map) != str:
                plt.scatter(combined_data['teff'], combined_data['logg'], 
                            label='Unresolved',
                            c=error_map, 
                            cmap='coolwarm', 
                            alpha=0.5, 
                            s=15)
                # Colorbar
                plt.colorbar()
            else:
                column = error_map
                plt.scatter(combined_data['teff'], combined_data['logg'], 
                            label='Unresolved',
                            c=combined_data[column], 
                            cmap='coolwarm', 
                            alpha=0.5, 
                            s=15)
                # Colorbar
                plt.colorbar(label=column)
        else:
            # Plot a scatter showing the unresloved and color by rchi2
            plt.scatter(combined_data['teff'], combined_data['logg'], 
                    label='Unresolved',
                    c=combined_data['rchi2'] * 1e3, 
                    cmap='coolwarm', 
                    alpha=0.5, 
                    s=15)  # Ensure a minimum size of 10
            # Colorbar
            plt.colorbar(label='rchi2')


    # plt.rcParams.update({'font.size': 20})
    plt.xlabel('teff')
    plt.ylabel('logg')
    plt.ylim(-1, 5)
    plt.legend()
    # Increase point size in legend
    legend = plt.gca().get_legend()
    for text in legend.get_texts():
        text.set_fontsize(20)
    for line in legend.get_lines():
        line.set_linewidth(5)
    plt.gca().invert_xaxis()
    plt.gca().invert_yaxis()
    plt.xlim(8000, 3000)
    plt.tight_layout()

    plt.show()


import seaborn as sns

def plot_hr_GALAH_heatmap(id=None, binary_data=None, lines=False, error_map=False, components=False):

    global GALAH_DR4
    global isochrone_data

    if GALAH_DR4 is None:
        print("Loading GALAH DR4 data into cache...")
        GALAH_DR4_dir = '/avatar/buder/GALAH_DR4/'
        GALAH_DR4 = df.FitsToDF(GALAH_DR4_dir + "catalogs/galah_dr4_allspec_240207.fits")

    if isochrone_data is None:
        print("Loading isochrone data into cache...")
        isochrone_data = Table.read('/avatar/yanilach/PhD-Home/binaries_galah-main/spectrum_analysis/BinaryAnalysis/assets/parsec_isochrones_reduced.fits')

    GALAH_data = GALAH_DR4

    plt.figure(figsize=(15, 7))
    cleaned_data = GALAH_data[['teff', 'logg']].replace([np.inf, -np.inf], np.nan).dropna()

    # Plot a smoothed heatmap instead of individual points
    sns.kdeplot(x=cleaned_data['teff'], y=cleaned_data['logg'], 
                cmap='magma', fill=True, levels=100, bw_adjust=0.5)

    age = 9.3
    isochrone_data_selection = isochrone_data[
        (isochrone_data['logAge'] <= age + 0.01) &
        (isochrone_data['logAge'] >= age - 0.01) &
        (isochrone_data['m_h'] == 0.0)
    ]

    ms_binaries = [(3000, 4.4), (3500, 4.4), (4000, 4.4), (4500, 4.4), 
                   (5000, 4.3), (5500, 4.2), (5750, 4)]
    min_teff = np.min([x for x, y in ms_binaries])
    max_teff = np.max([x for x, y in ms_binaries])

    # Fit and plot a polynomial to binary main sequence
    x, y = zip(*ms_binaries)
    coefficients = np.polyfit(x, y, 2)
    polynomial = np.poly1d(coefficients)
    plt.plot(np.arange(3500, max_teff), polynomial(np.arange(3500, max_teff)), 
             color='red', label='Binary Main Sequence', ls='--')

    plt.plot(10 ** isochrone_data_selection['logT'], 
             isochrone_data_selection['logg'], 
             alpha=0.5, color='black', label='MS', ls='dotted')

    plt.rcParams.update({'font.size': 20})
    plt.xlabel('Teff (K)')
    plt.ylabel('Log g')
    plt.ylim(-1, 5)
    plt.legend()
    plt.gca().invert_xaxis()
    plt.gca().invert_yaxis()
    plt.xlim(8000, 3000)
    plt.tight_layout()

    plt.show()
