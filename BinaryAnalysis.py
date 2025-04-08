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

import pyswarms as ps
from pyswarms.single.global_best import GlobalBestPSO
from pyswarms.utils.plotters import plot_cost_history
from pyswarms.utils.functions import single_obj as fx

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

from scipy.interpolate import LinearNDInterpolator
isochrone_table = Table.read('/avatar/yanilach/PhD-Home/binaries_galah-main/spectrum_analysis/BinaryAnalysis/' + 'assets/' + 'parsec_isochrones_extended.fits')
isochrone_interpolator = af.load_isochrones(type='trilinear')

sobject_id = int(sys.argv[1])
tmass_id = str(sys.argv[2])
labels_interp =  ['f_contr', 'mass_1', 'rv_1', 'vmic_1', 'vsini_1', 'mass_2', 'rv_2', 'vmic_2', 'vsini_2', 'FeH', 'age', 'teff_1', 'teff_2', 'logg_1', 'logg_2', 'logl_1', 'logl_2']


def fit_model(sobject_id):
    
    start_time = time.time()

    spectrum = af.read_spectrum(sobject_id, tmass_id)
    if spectrum == False:
        return

    # same_fe_h = False
    file = '/avatar/buder/GALAH_DR4/analysis_products_single/'+str(sobject_id)[:6]+'/'+str(sobject_id)+'/'+str(sobject_id)+'_single_fit_results.fits'

    try:
        single_results = Table.read(file)
    except:
        print('Single results not available')
        return
    

    # model = StellarModel(labels = ['teff', 'logg', 'rv', 'fe_h', 'vmic', 'vsini']) # Model with no interpolation
    # model = StellarModel(id=sobject_id, labels = ['mass', 'age', 'metallicity', 'rv', 'fe_h', 'vmic', 'vsini'], interpolator=isochrone_interpolator, interpolate_flux=True) # Flux can be used as a free parameter (False) or can be determined from luminosity ratios (from the isochrone) (True)
    
    
    # Use same FeH and age for both components. Use interpolation.
    # model = StellarModel(
    #             id=sobject_id, 
    #             labels = ['mass', 'rv', 'vmic', 'vsini'], 
    #             # fixed_labels=['FeH'], # These are parameters that exist, but we don't fit for. They are used as is or are outputs of the interpolation?
    #             single_labels=['FeH', 'age'], # These are parameters that are the same for both components
    #             interpolator=isochrone_interpolator, 
    #             interpolate_flux=True,
    #             same_fe_h=True
    #         ) # Flux can be used as a free parameter (False) or can be determined from luminosity ratios (from the isochrone) (True)


    # Use a model where we don't interpolate
    # model = StellarModel(
    #             id=sobject_id, 
    #             labels = ['teff', 'logg', 'mass', 'rv', 'vmic', 'vsini'], 
    #             # fixed_labels=['FeH'], # These are parameters that exist, but we don't fit for. They are used as is or are outputs of the interpolation?
    #             single_labels=['FeH', 'age'], # These are parameters that are the same for both components
    #             interpolate_flux=False,
    #             same_fe_h=True
    #         )

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

    # Trilinear Interpolation
    model = StellarModel(
                id=sobject_id, 
                labels = ['mass', 'rv', 'vmic', 'vsini'], 
                derived_labels=['teff', 'logg', 'logl'], # These are parameters that exist, but we don't fit for. They are used as is or are outputs of the interpolation?
                single_labels=['FeH', 'age'], # These are parameters that are the same for both components
                interpolator='trilinear', 
                isochrone_table=isochrone_table,
                interpolate_flux=True,
                same_fe_h=True
            ) # Flux can be used as a free parameter (False) or can be determined from luminosity ratios (from the isochrone) (True)

    model.bounds['f_contr'] = (0, 1)

    # Same bounds for both components. Overwrite with model.bounds['rv_1'] == x if required
    model.set_bounds('teff', (3, 8))
    model.set_bounds('logg', (0.0, 5.0))
    model.set_bounds('vmic', (0, 4))
    model.set_bounds('vsini', (0, 30))

    age_min = (10**isochrone_table['logAge'].min()) / 1e9
    age_max = (10**isochrone_table['logAge'].max()) / 1e9
    model.set_bounds('age', (age_min, age_max))
    

    model.set_bounds('FeH', (isochrone_table['m_h'].min(), isochrone_table['m_h'].max()))

    try:
        model.set_bounds('mass', (isochrone_table['mass'].min(), isochrone_table['mass'].max()))
    except:
        model.set_bounds('mass', (isochrone_table['mini'].min(), isochrone_table['mini'].max()))


    # model.set_bounds('metallicity', (isochrone_table['m_h'].min(), isochrone_table['m_h'].max()))
    model.set_bounds('logL', (isochrone_table['logL'].min(), isochrone_table['logL'].max()))

    model.params['f_contr'] = 0.5


    model.params['rv_1'] = single_results['rv_gauss'][0]
    
    # Assume the RVs are within 10 km/s of the values found by the CCF
    model.bounds['rv_1'] = (model.params['rv_1'] - 30, model.params['rv_1'] + 30)


    model.params['rv_2'] = single_results['rv_peak_2'][0]
    if np.isnan(model.params['rv_2']) or str(model.params['rv_2']) == "--":

        CCF_results = pd.read_csv("/avatar/yanilach/PhD-Home/binaries_galah-main/spectrum_analysis" + "/CCF_results_processed.txt", sep='\t', dtype={'no_peaks': float, 'rv_1': float, 'rv_2': float})
        CCF_results = CCF_results[CCF_results['sobject_id'] == sobject_id]
        # print("CCF results: ", CCF_results, CCF_results['rv_1'].values[0], CCF_results['rv_2'].values[0])
        ccf = True

        # If CCF_results matched to an object ID and there are more than 1 peaks, use the RV values from the CCF
        if len(CCF_results) > 0 and not np.isnan(CCF_results['rv_1'].values[0]) and not np.isnan(CCF_results['rv_2'].values[0]) and CCF_results['no_peaks'].values[0] > 1:
        # if CCF_results['no_peaks'].values[0] > 1:
            model.params['rv_1'] = CCF_results['rv_1'].values[0]
            model.params['rv_2'] = CCF_results['rv_2'].values[0]
            
            # Update both bounds as we've potentially changed both RV1 and RV2
            model.bounds['rv_1'] = (model.params['rv_1'] - 10, model.params['rv_1'] + 10)
            model.bounds['rv_2'] = (model.params['rv_2'] - 10, model.params['rv_2'] + 10)
            print("Using CCF results for RVs: ", model.params['rv_1'], model.params['rv_2'])

        # Use the RV1 value and assume RV2 is 100 km/s away. Best guess based on existing results.
        elif model.params['rv_1'] > 0:
            model.params['rv_2'] = model.params['rv_1'] - 100
            ccf = False
        else:
            model.params['rv_2'] = model.params['rv_1'] + 100
            ccf = False

        # Update the bounds
        if not ccf:
            model.bounds['rv_2'] = (model.params['rv_1'] - 200, model.params['rv_1'] + 200)
            print("No RV2 value! Setting rv_2 to ", model.params['rv_2'])

    else:
        model.bounds['rv_2'] = (model.params['rv_2'] - 30, model.params['rv_2'] + 30)

    model.set_param('teff', single_results['teff'][0]/1000.)
    model.set_param('logg', single_results['logg'][0])

    # Set the parameters required for the interpolation. The resulting interpolated parameters (teff, logg, and logL) are highly sensitive to these.
    model.set_param('age', float(sys.argv[3]))
    # Update the age bounds to be within +- 2 and constrained by the bounds of the isochrone
    model.set_bounds('age', (max(age_min, model.params['age'] - 2), min(age_max, model.params['age'] + 2)))

    mass_init = float(sys.argv[4])
    model.set_param('mass', mass_init)
    # Update the mass bounds to be within +- 2 and constrained by the bounds of the isochrone
    # model.set_bounds('mass', (max(isochrone_table['mass'].min(), mass_init - 2), min(isochrone_table['mass'].max(), mass_init + 2)))
    model.set_bounds('mass', (max(isochrone_table['mini'].min(), mass_init - 2), min(isochrone_table['mini'].max(), mass_init + 2)))

    model.set_param('FeH', float(sys.argv[5]))
    # Update FeH bounds to be within +- 1 and constrained by the isochrone
    # model.set_bounds('FeH', (max(isochrone_table['m_h'].min(), model.params['FeH'] - 1), min(isochrone_table['m_h'].max(), model.params['FeH'] + 1)))
    model.set_bounds('FeH', (max(isochrone_table['m_h'].min(), model.params['FeH'] - 1), min(isochrone_table['m_h'].max(), model.params['FeH'] + 1)))


    model.set_param('vmic', 1.5)
    model.set_param('vsini', 4.0)

    # Clip all the parameters to be within the bounds of the model
    for param in model.get_labels():
        # print("Clipping parameter ", param, model.params[param], model.bounds[param])
        # if param is not finite or param is None set to middle of the bounds
        if not np.isfinite(model.params[param]) or model.params[param] is None:
            model.params[param] = 1
        
        model.params[param] = np.clip(model.params[param], model.bounds[param][0], model.bounds[param][1])

    af.load_neural_network(spectrum)
    af.set_iterations(0)
    af.load_dr3_lines()
    
    # Generate an initial model with the starting parameters
    print("Initial parameters:")
    # print(model.params)



    print("Initial bounds:")
    print(model.bounds)
    model.unnormalized_bounds = model.bounds

    print("Unormalized bounds:")
    print(model.unnormalized_bounds)
    # This is not just for the output: it sets the parameters in the model object for the intialisation. DO NOT REMOVE.
    model.interpolate()
    print(model.get_params(values_only=False, exclude_fixed=True))
    # print(model.fixed_labels)

    if len(model.get_params(values_only=True, exclude_fixed=True)) != len(model.bounds):
        print("Length of parameters and bounds do not match")
        # Print out which parameters don't have bounds
        for p in model.get_params(values_only=False, exclude_fixed=True):
            if p not in model.bounds:
                print("Missing bound: ", p)
        print("Parameters: ", model.get_params(values_only=False, exclude_fixed=True))
        print("Length of parameters and bounds mismatch", len(model.get_params(values_only=True, exclude_fixed=True)), len(model.bounds))
        return
    else:
        print("Length of parameters and bounds match", len(model.get_params(values_only=True, exclude_fixed=True)), len(model.bounds))

    wave_init, data_init, sigma2_init, model_init, unmasked_init = af.return_wave_data_sigma_model(model, spectrum, model.same_fe_h)
    unmasked = unmasked_init

    # Produce a plot with the initial parameters
    model.generate_model(spectrum)
    model.plot()



    iterations = 0
    def objective_function_norm(normalized_params):

        global previous_params
        global iterations
        show_plot = True
        residuals_list = []

        # If this is LBFGS-B, there is only a single particle. If PSO we need to iterate over all particles
        if isinstance(normalized_params[0], (list, np.ndarray)):
            show_plot = False
            normalized_params_list = normalized_params
        else:
            show_plot = True
            normalized_params_list = [normalized_params]

        # Note we pass in normalized parameters. We need to denormalize them before passing them to the model
        for i, particle_params in enumerate(normalized_params_list):

            # Denormalize the parameters
            model_parameters = denormalize_parameters(particle_params, original_bounds)

            # Calculate the model flux using the current parameters
            model_flux = af.get_flux_only(wave_init, model, spectrum, model.same_fe_h, unmasked, *model_parameters, plot=False)
            
            # Need to generate a model with the current parameters to determine residual
            model.generate_model(spectrum)
            # residuals = model.get_residual()

            in_bounds = True
            # for param, bounds in zip( model.get_params(), unnormalized_bounds):
            #     # print("Checking bounds for ", param, model.params[param], bounds)
            #     if model.params[param] < bounds[0] or model.params[param] > bounds[1]:
            #         # print(f"Parameter {param} out of bounds: {model.params[param]}. Bounds: {bounds[0]} and {bounds[1]}")
            #         residuals = 1e10
            #         residuals_list.append(residuals)
            #         in_bounds = False
            #         break

            if in_bounds:
                residuals = model.get_rchi2()
                residuals_list.append(residuals)

            previous_params = np.copy(particle_params)

        best_params = optimizer.swarm.best_pos
        # if len(best_params) > 0:
        #     best_params = denormalize_parameters(best_params, original_bounds)
        #     model.set_params(best_params)
        #     print(f"Global Best Parameters: {best_params}")
        #     print(model.params)
        #     model.generate_model(spectrum)
        #     model.plot(title_text=str(iterations))
            # iterations += 1

        # return residuals
        if len(residuals_list) == 1:
            return residuals_list[0]
        else:
            return residuals_list

# 
    try:
        # Fit the model to the data. This takes the model parameters and produces a synthetic spectra using the neural network. It then compares this to the observed data and adjusts the model parameters (and thereby the synthetic spectra from the NN) to minimize the difference between the two.
        kwargs={'maxfev':20,'xtol':1e-5, 'gtol':1e-5, 'ftol':1e-5}
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
    except Exception as e:
        print(f"Curve fitting failed with error: {e}")
        pass




    params_curve_fit = model.get_params(values_only=True)
    params_curve_fit_list = ', '.join(map(str, params_curve_fit))
    curve_fit_rchi2 = model.get_rchi2(with_bounds=False)

    # Get the original parameter values and bounds
    original_params = model.get_params(values_only=True)# model_parameters_iter1 # model.get_params(values_only=True)
    original_bounds = model.get_bounds(type='tuple')
    model.original_bounds = original_bounds
    bounds = model.get_bounds(type='tuple')


    # Normalize the initial parameter values
    normalized_x0 = normalize_parameters(original_params, bounds)
    
    # If the curve fit hasn't failed, assume the parameters returned from it are good. 
    # Use these new values to update the bounds for the minimization and as the initial guess for the next optimization.
    if model.get_residual() < 5:
        # print("Good curve fit. Updating bounds and initial parameters")
        # Update bounds to so that they are within 25% of the values returned by curve_fit
        margin = 0.2
            
    else:
        # Normalize the bounds
        # bounds = [(0, 1)] * len(bounds)
        # unnormalized_bounds = original_bounds
        margin = 1


    bounds = [(max(lb, p - margin * abs(p)), min(ub, p + margin * abs(p)) ) for p, (lb, ub) in zip(normalized_x0, [(0, 1)] * len(bounds))]

    # Set a 0.1 margin on the bounds required for the interpolation. Otherwise we create an interpolator that is too large if the margin on the parameters is to expansive.
    if margin == 1:
        bounds[model.label('age')] = (max(bounds[model.label('age')][0], normalized_x0[model.label('age')] - 0.1), min(bounds[model.label('age')][1], normalized_x0[model.label('age')] + 0.1))
        bounds[model.label('FeH')] = (max(bounds[model.label('FeH')][0], normalized_x0[model.label('FeH')] - 0.1), min(bounds[model.label('FeH')][1], normalized_x0[model.label('FeH')] + 0.1))
        bounds[model.label('mass', comp=1)] = (max(bounds[model.label('mass', comp=1)][0], normalized_x0[model.label('mass', comp=1)] - 0.1), min(bounds[model.label('mass', comp=1)][1], normalized_x0[model.label('mass', comp=1)] + 0.1))


    # Override these bounds
    bounds[model.label('f_contr')] = [0, 1]
    bounds[model.label('teff', 1)] = [0, 1]
    bounds[model.label('teff', 2)] = [0, 1]
    bounds[model.label('logg', 1)] = [0, 1]
    bounds[model.label('logg', 2)] = [0, 1]
    # Ensure these new bounds are also within the (0, 1) bounds
    bounds = [(max(0, lb), min(1, ub)) for lb, ub in bounds]




    unnormalized_bounds = [(original_lb + (lb - 0) * (original_ub - original_lb), 
                            original_lb + (ub - 0) * (original_ub - original_lb)) 
                        for (lb, ub), (original_lb, original_ub) in zip(bounds, original_bounds)]
    print("Age bounds")
    print(unnormalized_bounds[model.label('age')])
    print("FeH bounds")
    print(unnormalized_bounds[model.label('FeH')])
    print("Mass bounds")
    print(unnormalized_bounds[model.label('mass', comp=1)])
    




    # Set the bounds of the model object. THIS IS IMPORTANT, we need to pass normalised bounds AS WELL as values.
    bounds_dict = {param: bound for param, bound in zip(model.get_params(values_only=False, exclude_fixed=True), bounds)}
    model.bounds = bounds_dict
    
    unnormalized_bounds_dict = {param: bound for param, bound in zip(model.get_params(values_only=False, exclude_fixed=True), unnormalized_bounds)}
    model.unnormalized_bounds = unnormalized_bounds_dict

    # print(unnormalized_bounds, model.label('age'))
    age_range = np.log10(np.array(unnormalized_bounds[model.label('age')]) * 1e9)
    m_h_range = unnormalized_bounds[model.label('FeH')]
    mass_range = unnormalized_bounds[model.label('mass', comp=1)]

    # Try creating the interpolator with the original range
    max_expansions = 10
    expansion_count = 0
    while expansion_count <= max_expansions:
        try:
            interpolator = af.create_interpolator(isochrone_table, age_range, m_h_range, mass_range)
            print(f"Interpolator created successfully after {expansion_count} expansions.")
            break
        except Exception as e:
            print(f"Interpolator creation failed on attempt {expansion_count + 1}: {e}")
            # Expand age_range, m_h_range, and mass_range
            if age_range[1] - age_range[0] < 1:
                delta = 1 + expansion_count * 0.1
                age_range = np.array([model.params['age'] - delta, model.params['age'] + delta])
                age_range = np.log10(age_range * 1e9)
                # Ensure the age range is within the bounds of the isochrone table
                age_range = np.clip(age_range, isochrone_table['logAge'].min(), isochrone_table['logAge'].max())
                print(f"Expanding age range bounds to attempt again. New age range: {age_range}")
            else:
                print("Age range is already sufficiently large, attempting again without further expansion.")

            delta = 0.1 + expansion_count * 0.1
            if m_h_range[1] - m_h_range[0] < 0.5:
                m_h_range = [m_h_range[0] - delta, m_h_range[1] + delta]
                # Ensure the m_h range is within the bounds of the isochrone table
                m_h_range = [max(m_h_range[0], isochrone_table['m_h'].min()), min(m_h_range[1], isochrone_table['m_h'].max())]
                print(f"Expanding m_h range bounds to attempt again. New m_h range: {m_h_range}")
            else:
                print("m_h range is already sufficiently large, attempting again without further expansion.")

            if mass_range[1] - mass_range[0] < 0.5:
                mass_range = [mass_range[0] - delta, mass_range[1] + delta]
                # Ensure the mass range is within the bounds of the isochrone table
                mass_range = [max(mass_range[0], isochrone_table['mini'].min()), min(mass_range[1], isochrone_table['mini'].max())]
                print(f"Expanding mass range bounds to attempt again. New mass range: {mass_range}")
            else:
                print("Mass range is already sufficiently large, attempting again without further expansion.")

            expansion_count += 1

    if expansion_count > max_expansions:
        raise ValueError(f"Failed to create an interpolator after {max_expansions} expansions.")

    print("Finished creating interpolator")


    model.interpolator = interpolator


    print("Optimizing with PSO")
    print("Params:")
    print(model.get_params())
    print("Updated bounds: ")
    print(unnormalized_bounds)
    print("Normalised params:")
    print(normalized_x0)
    print("Normalised bounds: ")
    print(model.get_bounds(type='tuple'))

    # result = scipy.optimize.minimize(
    #     objective_function_norm,
    #     x0=normalized_x0, #model.get_params(values_only=True),
    #     method='L-BFGS-B',
    #     bounds= bounds, #[(0, 1)] * len(bounds), #model.get_bounds(type='tuple'),
    #     # Ftol is the relative error desired in the sum of squares.
    #     # Gtol is the gradient norm desired in the sum of squares.
    #     options={'maxfun': 10000, 'gtol': 1e-10, 'ftol': 1e-10, 'eps': 1e-5}
    # )

    # lower_bounds = [b[0] for b in bounds]
    # upper_bounds = [b[1] for b in bounds]

    # # Call the pso function
    # best_params, best_score = pso(
    #     objective_function_norm,
    #     lower_bounds,
    #     upper_bounds,
    #     swarmsize=100,        # Number of particles
    #     maxiter=50,          # Number of iterations
    #     minstep=1e-8,         # Minimum step size before convergence
    #     minfunc=1e-8,         # Minimum function change before convergence
    #     debug=True            # Set to True for detailed logging
    # )

    lower_bounds = [b[0] for b in bounds]
    upper_bounds = [b[1] for b in bounds]
    bounds_array = [lower_bounds, upper_bounds]

    # print("Bounds array", bounds_array)
    # return

    # Define options
    options = {
        'c1': 0.8,  # Cognitive parameter
        'c2': 0.3,  # Social parameter
        'w': 0.7,   # Inertia weight OR 0.8
    }

    # Initialize the optimizer
    optimizer = ps.single.GlobalBestPSO(
        n_particles=100,              # Number of particles
        dimensions=len(lower_bounds),       # Dimensionality of the parameter space
        options=options,
        bounds=bounds_array           # (lower_bounds, upper_bounds)
    )

    # Perform optimization
    best_score, best_params = optimizer.optimize(
        objective_function_norm,     # Objective function
        iters=60,                   # Number of iterations
        verbose=True,                 # Display progress    
    )

    # correct this! original_bounds
    final_params = denormalize_parameters(best_params, original_bounds)

    for i, param in enumerate(model.get_params(values_only=False)):
        model.params[param] = final_params[i]
    
    # Regenerate the model with the final parameters
    model.interpolate()
    model.generate_model(spectrum)

    params = model.get_params(values_only=True)
    params_list = ', '.join(map(str, params))

    total_time = time.time() - start_time

    interpolated_params = model.interpolate(return_values=True)
    interpolated_params = ', '.join(map(str, interpolated_params[0])) + ', ' + ', '.join(map(str, interpolated_params[1]))

    print(model.get_residual(), model.get_rchi2(with_bounds=False), params_list, curve_fit_rchi2, params_curve_fit_list, interpolated_params, total_time)

fit_model(sobject_id)