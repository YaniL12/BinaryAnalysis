# Basic packages
import numpy as np
import time
import sys
import os
from pathlib import Path
import logging
import importlib
# import mysql.connector

import pandas as pd
pd.set_option('display.max_columns', None)

# Astropy packages
from astropy.table import Table
from astropy.io import fits

# Matplotlib packages
import matplotlib.pyplot as plt

# Scipy
import scipy
from scipy.optimize import curve_fit
from scipy import signal

import pyswarm
from pyswarm import pso


import warnings
warnings.filterwarnings('ignore', category=UserWarning, append=True)
warnings.filterwarnings('ignore')


# Useful if working in SSH Vscode
working_directory = '/avatar/yanilach/PhD-Home/binaries_galah-main/spectrum_analysis/BinaryAnalysis'
os.chdir(working_directory)
import AnalysisFunctions as af
from stellarmodel import StellarModel

def normalize_parameters(params, bounds):
    normalized_params = [(p - lb) / (ub - lb) for p, (lb, ub) in zip(params, bounds)]
    return normalized_params

def denormalize_parameters(normalized_params, bounds):
    denormalized_params = [lb + n * (ub - lb) for n, (lb, ub) in zip(normalized_params, bounds)]
    return denormalized_params



sys.path.append(os.path.join(working_directory, 'utils'))
import AstroPandas as ap

isochrone_table = Table.read(working_directory +  '/assets/parsec_isochrones_logt_8p00_0p01_10p17_mh_m2p75_0p25_m0p75_mh_m0p60_0p10_0p70_GaiaEDR3_2MASS.fits')
isochrone_interpolator = af.load_isochrones()

sobject_id = int(sys.argv[1])
tmass_id = str(sys.argv[2])



def fit_model(sobject_id):
    spectrum = af.read_spectrum(sobject_id, tmass_id)
    if spectrum == False:
        return

    # same_fe_h = False

    try:
        single_results = Table.read('/avatar/buder/GALAH_DR4/analysis_products_single/'+str(sobject_id)[:6]+'/'+str(sobject_id)+'/'+str(sobject_id)+'_single_fit_results.fits')
    except:
        print('Single results not available')
        return

    # model = StellarModel(labels = ['teff', 'logg', 'rv', 'fe_h', 'vmic', 'vsini']) # Model with no interpolation
    # model = StellarModel(id=sobject_id, labels = ['mass', 'age', 'metallicity', 'rv', 'fe_h', 'vmic', 'vsini'], interpolator=isochrone_interpolator, interpolate_flux=True) # Flux can be used as a free parameter (False) or can be determined from luminosity ratios (from the isochrone) (True)
    
    
    # Use same FeH and age for both components
    model = StellarModel(
                id=sobject_id, 
                labels = ['mass', 'rv', 'vmic', 'vsini'], 
                # fixed_labels=['FeH'], # These are parameters that exist, but we don't fit for. They are used as is or are outputs of the interpolation?
                single_labels=['FeH', 'age'], # These are parameters that are the same for both components
                interpolator=isochrone_interpolator, 
                interpolate_flux=True,
                same_fe_h=True
            ) # Flux can be used as a free parameter (False) or can be determined from luminosity ratios (from the isochrone) (True)

    # # Use different FeH for each component and allow age to vary
    # model = StellarModel(
    #             id=sobject_id, 
    #             labels = ['mass', 'rv', 'vmic', 'vsini', 'FeH', 'age'], 
    #             # fixed_labels=['FeH'], # These are parameters that exist, but we don't fit for. They are used as is or are outputs of the interpolation?
    #             # single_labels=['FeH', 'age'], # These are parameters that are the same for both components
    #             interpolator=isochrone_interpolator, 
    #             interpolate_flux=True,
    #             same_fe_h=False
    #         ) # Flux can be used as a free parameter (False) or can be determined from luminosity ratios (from the isochrone) (True)

    model.bounds['f_contr'] = (0, 1)

    # Same bounds for both components. Overwrite with model.bounds['rv_1'] == x if required
    model.set_bounds('teff', (3, 8))
    model.set_bounds('logg', (0.0, 5.0))
    model.set_bounds('vmic', (0, 4))
    model.set_bounds('vsini', (0, 30))

    age_min = (10**isochrone_table['logAge'].min()) / 1e9
    age_max = (10**isochrone_table['logAge'].max()) / 1e9

    model.set_bounds('FeH', (-4.0, 1.0))
    model.set_bounds('age', (age_min, age_max))
    model.set_bounds('mass', (isochrone_table['mass'].min(), isochrone_table['mass'].max()))
    # model.set_bounds('metallicity', (isochrone_table['m_h'].min(), isochrone_table['m_h'].max()))
    model.set_bounds('logL', (isochrone_table['logL'].min(), isochrone_table['logL'].max()))

    model.params['f_contr'] = 0.5


    model.params['rv_1'] = single_results['rv_gauss'][0]
    model.params['rv_2'] = single_results['rv_peak_2'][0]
    if np.isnan(model.params['rv_2']) or str(model.params['rv_2']) == "--":
        print("No RV2 value!")
        return

    min_rv = min(model.params['rv_1'], model.params['rv_2']) - 100
    max_rv = max(model.params['rv_1'], model.params['rv_2']) + 100
    model.set_bounds('rv', (min_rv, max_rv))

    model.set_param('teff', single_results['teff'][0]/1000.)
    model.set_param('logg', single_results['logg'][0])

    # Set the parameters required for the interpolation
    model.set_param('age', float(sys.argv[3]))
    model.set_param('mass', float(sys.argv[4]))
    model.set_param('FeH', float(sys.argv[5]))

    model.set_param('vmic', 1.5)
    model.set_param('vsini', 4.0)

    af.load_neural_network(spectrum)
    af.set_iterations(0)
    af.load_dr3_lines()
    
    # Generate an initial model with the starting parameters
    print("Initial parameters:")
    # print(model.params)
    print(model.get_params(values_only=False, exclude_fixed=True))
    # print(model.fixed_labels)

    print("Initial bounds:")
    print(model.bounds)


    if len(model.get_params(values_only=True, exclude_fixed=True)) != len(model.bounds):
        print("Length of parameters and bounds do not match")
        # Print out which parameters don't have bounds
        for p in model.get_params(values_only=False, exclude_fixed=True):
            if p not in model.bounds:
                print("Missing: ", p)
        return
    else:
        print("Length of parameters and bounds match", len(model.get_params(values_only=True, exclude_fixed=True)), len(model.bounds))

    wave_init, data_init, sigma2_init, model_init, unmasked_init = af.return_wave_data_sigma_model(model, spectrum, model.same_fe_h)
    unmasked = unmasked_init

    # Produce a plot with the initial parameters
    model.generate_model(spectrum)
    model.plot()



    def objective_function_norm(normalized_params):

        global previous_params

        # Denormalize the parameters
        model_parameters = denormalize_parameters(normalized_params, model.get_bounds(type='tuple'))

        # Calculate the model flux using the current parameters
        model_flux = af.get_flux_only(wave_init, model, spectrum, model.same_fe_h, unmasked, *model_parameters, plot=False)
        
        # Need to generate a model with the current parameters to determine residual
        model.generate_model(spectrum)
        # residuals = model.get_residual()
        residuals = model.get_rchi2()


        # print('Step ', np.array(normalized_params - previous_params))
        previous_params = np.copy(normalized_params)
        
        return residuals

    # Fit the model to the data. This takes the model parameters and produces a synthetic spectra using the neural network. It then compares this to the observed data and adjusts the model parameters (and thereby the synthetic spectra from the NN) to minimize the difference between the two.
    kwargs={'maxfev':20000,'xtol':1e-5, 'gtol':1e-5, 'ftol':1e-5}
    model_parameters_iter1, covariances_iter1 = curve_fit(
        lambda wave_init, 
            *model_parameters: af.get_flux_only(wave_init, model, spectrum, model.same_fe_h, unmasked, *model_parameters, plot=True),
        wave_init[unmasked_init],
        data_init[unmasked_init],
        p0=model.get_params(values_only=True, exclude_fixed=True),
        sigma=np.sqrt(sigma2_init[unmasked_init]),
        absolute_sigma=True,
        bounds=model.get_bounds(),
        **kwargs
    )




    # Get the original parameter values and bounds
    original_params = model.get_params(values_only=True)# model_parameters_iter1 # model.get_params(values_only=True)
    bounds = model.get_bounds(type='tuple')


    # Normalize the initial parameter values
    normalized_x0 = normalize_parameters(original_params, bounds)
    
    # If the curve fit hasn't failed, assume the parameters returned from it are good. 
    # Use these new values to update the bounds for the minimization and as the initial guess for the next optimization.
    if model.get_residual() < 10:
        # print("Good curve fit. Updating bounds and initial parameters")
        # Update bounds to so that they are within 25% of the values returned by curve_fit
        margin = 0.1
        bounds = [(max(lb, p - margin * abs(p)), min(ub, p + margin * abs(p)) ) for p, (lb, ub) in zip(normalized_x0, [(0, 1)] * len(bounds))]
        # Ensure these new bounds are also within the (0, 1) bounds
        bounds = [(max(0, lb), min(1, ub)) for lb, ub in bounds]
            
    else:
        bounds = [(0, 1)] * len(bounds)


    # result = scipy.optimize.minimize(
    #     objective_function_norm,
    #     x0=normalized_x0, #model.get_params(values_only=True),
    #     method='L-BFGS-B',
    #     bounds= bounds, #[(0, 1)] * len(bounds), #model.get_bounds(type='tuple'),
    #     # Ftol is the relative error desired in the sum of squares.
    #     # Gtol is the gradient norm desired in the sum of squares.
    #     options={'maxfun': 10000, 'gtol': 1e-10, 'ftol': 1e-10, 'eps': 1e-5}
    # )

    lower_bounds = [b[0] for b in bounds]
    upper_bounds = [b[1] for b in bounds]

    # Call the pso function
    best_params, best_score = pso(
        objective_function_norm,
        lower_bounds,
        upper_bounds,
        swarmsize=100,        # Number of particles
        maxiter=200,          # Number of iterations
        minstep=1e-8,         # Minimum step size before convergence
        minfunc=1e-8,         # Minimum function change before convergence
        debug=True            # Set to True for detailed logging
    )

    final_params = denormalize_parameters(best_params, model.get_bounds(type='tuple'))

    for i, param in enumerate(model.get_params(values_only=False)):
        model.params[param] = final_params[i]
        
    params = model.get_params(values_only=True)
    params_list = ', '.join(map(str, params))
    print(model.get_residual(), model.get_rchi2(), params_list)


fit_model(sobject_id)