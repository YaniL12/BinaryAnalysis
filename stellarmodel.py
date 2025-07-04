import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import AnalysisFunctions as af
from pandas import DataFrame
from pandas import Series as Series
import time
import os
important_lines, important_molecules = af.load_dr3_lines()
from astropy.table import Table

import sys
working_directory = '/avatar/yanilach/PhD-Home/binaries_galah-main/spectrum_analysis/BinaryAnalysis'
sys.path.append(os.path.join(working_directory, 'utils'))
import AstroPandas as ap

class StellarModel:
    interpolator = None

    # Constructor with default labels and components. Override for custom labels and components
    def __init__(self, id="No ID", labels = ['rv', 'teff', 'logg', 'FeH', 'vmic', 'vsini'], derived_labels=[], fixed_labels=[], single_labels=[], components = 2, same_fe_h=False, interpolator=None, interpolate_flux=False, original_bounds=None, isochrone_table=None, load_from_results=False):

    # Dictionaries for labels, bounds, and initial values
    # These should be instance level variables not class level. Otherwise they will be shared between all instances.
        self.id = id
        self.components = components
        self.same_fe_h = same_fe_h
        self.spectrum = None

        self.model_labels = {}
        self.derived_labels = derived_labels
        self.unique_labels = []
        self.single_labels = single_labels

        self.bounds = {}
        self.original_bounds = original_bounds
        self.unnormalized_bounds = {}
        self.params = {}
        self.indices = {}
        self.unique_indices = {}
        self.wavelengths = []
        self.flux = []
        self.flux_uncertainty = []
        self.model_flux = []
        self.component_model_fluxes = {}

        self.param_data = {}

        self.isochrone_table = isochrone_table

        # Should we use a simple flux model or determine flux ratios from interpolated luminosity?
        self.interpolate_flux = interpolate_flux
        if interpolate_flux and interpolator is None:
            raise ValueError("You have set interpolate_flux to True but have not provided an interpolator. Please provide an interpolator to use this feature.")

        self.fixed_labels = fixed_labels

        # self.unique_labels.append('f_contr')
        self.unique_labels.extend(labels)
        self.unique_labels.extend(derived_labels)

        # Only one instance of f_contr. This will need to be changeed for multiple systems where n_comps > 2
        self.model_labels['f_contr'] = 'f_contr'

        for i, label in enumerate(labels):
            self.unique_indices[label] = i

        for j in range(self.components):
            for i, label in enumerate(list(labels) + list(derived_labels)):
                self.model_labels[label + '_' + str(j+1)] = label + '_' + str(j+1)

        # Add single labels to the model labels. These are those that are the same for both components.
        for i, label in enumerate(single_labels):
            self.model_labels[label] = label
            
        # Instantiate bounds and params dictionaries with defaults
        for i, label in enumerate(list(self.model_labels)):
            self.bounds[label] = (-1e5, 1e5)

            if label in self.model_labels:
                self.params[label] = 0
                self.indices[label] = i


        if interpolator is not None:
            self.interpolator = interpolator

            # Interpolator requires mass, age, and metallicity to be present in the model labels
            # Mass will be different for each component, the remaining parameters will be the same for both components (age, FeH).
            # if same_fe_h and all(label in self.single_labels for label in ['age', 'FeH']) and all(label in self.unique_labels for label in ['mass']):
            #     self.add_param('teff', 0)
            #     self.add_param('logg', 0)
            #     self.add_param('logl', 0)

            # # elif all(label in self.unique_labels for label in ['mass', 'age']) and all(label in self.fixed_labels for label in ['FeH']):
            # elif all(label in self.unique_labels for label in ['mass', 'age', 'FeH']):
            #     self.add_param('teff', 0)
            #     self.add_param('logg', 0)
            #     self.add_param('logl', 0)

                # self.interpolate()
        
        self.param_data = {key: [] for key in self.params.keys()}
        self.param_data['residual'] = []


        galah_dr4_directory = '/avatar/buder/GALAH_DR4/'
        dir = galah_dr4_directory + 'analysis_products_allstar/' + str(self.id)[:6] + '/' + str(self.id) + '/' + str(self.id) + '_allstar_fit_spectrum.fits'
        dir2 = galah_dr4_directory + 'analysis_products_single/' + str(self.id)[:6] + '/' + str(self.id) + '/' + str(self.id) + '_single_fit_spectrum.fits'
        dir3 = galah_dr4_directory + 'analysis_products_binary/' + str(self.id)[:6] + '/' + str(self.id) + '/' + str(self.id) + '_binary_fit_spectrum.fits'

        try:
            if os.path.exists(dir):
                fits_file = ap.FitsToDFWithVariableLengthCols(dir)
                fits_file = fits_file[0]
            elif os.path.exists(dir2):
                fits_file = ap.FitsToDFWithVariableLengthCols(dir2)
                fits_file = fits_file[0]
                print("Succsefully found file for object at analysis_products_single: " + dir2)
            elif os.path.exists(dir3):
                fits_file = ap.FitsToDFWithVariableLengthCols(dir3)
                fits_file = fits_file[0]
                print("Succsefully found file for object at analysis_products_binary: " + dir3)
            else:
                print("File not found: ", dir)
                raise FileNotFoundError(f"File not found: {dir}")
        except:
            print("No file found for spectra ", dir)

        self.flux_uncertainty = fits_file['uob'].values
        self.model_flux_single = fits_file['smod'].values
        self.wavelengths_single = fits_file['wave'].values

        self.rchi2_single = np.median(np.abs(fits_file['sob'] - self.model_flux_single)/ self.flux_uncertainty)

        # Load in spectral data. e.g. wavelengths, flux, and model flux.
        self.spectrum = af.read_spectrum(self.id)
        # May need to remove this if loading in data manually from a result file.
        if not load_from_results:
            self.generate_model()

        # Flux uncertainty may be slightly different in length from the model flux by a few points.
        # Interpolate the flux uncertainty to match the model flux
        print("Lenght of wavelengths: ", len(self.wavelengths))
        if self.wavelengths is not None and len(self.wavelengths) > 0:
            print("Interpolating flux uncertainty to match model flux")
            self.flux_uncertainty = np.interp(self.wavelengths, self.wavelengths_single, self.flux_uncertainty)
            # Need this for the CCF plotting.
            self.model_flux_single_original = self.model_flux_single
            self.model_flux_single = np.interp(self.wavelengths, self.wavelengths_single, self.model_flux_single)

    def save_data(self):
        for i, param in enumerate(self.params):
            self.param_data[param].append(self.params[param])

        r = 100 * np.sum(abs(np.array(self.model_flux) - np.array(self.flux))) / len(self.flux)
        self.param_data['residual'].append(r)

    def load_data(self, data):
        # print(type(data))
        # print("Type: ", type())
        # Check if data is a DataFrame
        if type(data) is DataFrame:
            row = data[data['sobject_id'] == self.id]
            if row.empty:
                print("ID not found in data. Check the ID you passed matches the model ID")
                return
            for col in data.columns:
                # print("Checking col: ", col)
                # If col exists in the model.params, set the value
                if col in self.get_labels() or col in self.derived_labels:
                    # print("Setting param: ", col, " to ", row[col].iloc[0])
                    self.params[col] = row[col].iloc[0]

            #Note we requre this as we interpolate for M = m - 0.1. REMOVE IF WE CHANGE THE OUTPUT FROM THE ANALYSIS!
            self.params['mass_1'] += 0.1
            self.params['mass_2'] += 0.1

            # self.params['age'] = 10 ** self.params['age'] / 1e9

            # Check for missing parameters that have not been set
            for col in self.get_labels() or col in self.derived_labels:
                if col not in data.columns:
                    print("Missing parameter: ", col)
                    # Raise an error
                    raise ValueError("Missing parameter: " + col + ". Please check the data you passed to the model is the original data table, for example not a merged table where parameters have prefixes.")
            
            af.load_neural_network(self.spectrum)
            
            print("Generating model")
            self.generate_model()

        elif isinstance(data, Series):
            for col in data.keys():
                # If col exists in the model.params, set the value
                if col in self.get_labels():
                    # print(col, data[col])
                    self.params[col] = data[col]
            # print(self.params)
            af.load_neural_network(self.spectrum)
            self.generate_model()

        if self.wavelengths is not None and len(self.wavelengths) > 0:
            print("Interpolating flux uncertainty to match model flux")
            self.flux_uncertainty = np.interp(self.wavelengths, self.wavelengths_single, self.flux_uncertainty)
            # Need this for the CCF plotting.
            self.model_flux_single_original = self.model_flux_single
            self.model_flux_single = np.interp(self.wavelengths, self.wavelengths_single, self.model_flux_single)

    # Call this if we have modified one or more parameters manually. This will update the model flux.
    # Effectively a load data call, without loadiing paramaters.
    def regenerate_spectrum(self):
        spectrum = af.read_spectrum(self.id)
        self.generate_model()
        self.spectrum = spectrum

    # def generate_interpolator(self):
        # age_range = np.log10(np.array(self.unnormalized_bounds['age']) * 1e9)
        # m_h_range = self.unnormalized_bounds['FeH']
        # mass_range = self.unnormalized_bounds['mass_1']

        # print("ranges ", age_range, m_h_range, mass_range)
        # print(self.isochrone_table)

        # interpolator = af.create_interpolator(self.isochrone_table, age_range, m_h_range, mass_range)

        # self.interpolator = interpolator
        
        # return interpolator




    def interpolate(self, return_values=False):
        # Accepts mass, log(age), metallicity.
        # Outputs Teff, logg, and log(L) bolometric (flux)

        # start_time = time.time()

        # This interpolation function is used before the initial curve fit. This creates an interpolator and interpolates on the fly. 
        # Code is also stored within AnalysisFunctions.py.
        if self.interpolator == 'trilinear':
            # if return_values:
                # print("Using trilinear interpolation")

            age_query = np.log10(self.params['age'] * 1e9)
            m_h_query = self.params['FeH']
            mass_query_1 = self.params['mass_1'] - 0.1
            mass_query_2 = self.params['mass_2'] - 0.1

            # print("Passing", self.params['mass_1'], np.log10(self.params['age'] * 1e9), self.params['FeH'], "\n")
            # print("Passing to analysis functions: ", mass_query_1, age_query, m_h_query)

            interpolate_1 = af.interpolate_isochrone(mass_query_1, age_query, m_h_query)
            interpolate_2 = af.interpolate_isochrone(mass_query_2, age_query, m_h_query)

            # print("Recieved", interpolate_1, "\n")


            self.params['teff_1'] = (interpolate_1['teff']) / 1000
            self.params['logg_1'] = interpolate_1['logg']
            self.params['logl_1'] = interpolate_1['logl']

            self.params['teff_2'] = (interpolate_2['teff'])  / 1000
            self.params['logg_2'] = interpolate_2['logg']
            self.params['logl_2'] = interpolate_2['logl']


        # After the curve fit we perform all subsequent interpolations using the cached interpolator generated in the main script by AnalysisFunctions.
        # This has reduced bounds and should be faster.
        elif self.interpolator is not None:
            # print("Cached interpolation.")
            # if return_values:
                # print("Using internal interpolation")

            # if all(label in self.unique_labels for label in ['mass', 'age']) and all(label in self.fixed_labels for label in ['FeH']):
                # Both components will have the same starting values.
                # Provide the log of the age in Gyr for interpolation
                # TODO Change this to be dynamic based on the number of components

            if self.same_fe_h:
                # Star 1
                interpolate_1 = self.interpolator(self.params['mass_1'], np.log10(self.params['age'] * 1e9), self.params['FeH'])
                # Star 2
                interpolate_2 = self.interpolator(self.params['mass_2'], np.log10(self.params['age'] * 1e9), self.params['FeH'])
            else:
                # Star 1
                interpolate_1 = self.interpolator(self.params['mass_1'], np.log10(self.params['age_1'] * 1e9), self.params['FeH_1'])
                # Star 2
                interpolate_2 = self.interpolator(self.params['mass_2'], np.log10(self.params['age_2'] * 1e9), self.params['FeH_2'])

            # self.set_param('teff', (10 ** interpolate[0]) / 1000)
            # self.set_param('logg', interpolate[1])
            # self.set_param('logl', interpolate[2])

            self.params['teff_1'] = (10 ** interpolate_1[0]) / 1000
            self.params['logg_1'] = interpolate_1[1]
            self.params['logl_1'] = interpolate_1[2]

            self.params['teff_2'] = (10 ** interpolate_2[0]) / 1000
            self.params['logg_2'] = interpolate_2[1]
            self.params['logl_2'] = interpolate_2[2]

        # The optimiser has gone outside the bounds. Set parmaeters to unreasonable values. This should results in a high residual.
        # Consider scaling parameters to prevent this in the optimiser (TODO)
        if np.any(np.isnan(self.params['teff_1'])) or np.any(np.isnan(self.params['teff_2'])):
            # print("Interpolated values are NaN. Check input values for interpolation")
            self.params['teff_1'] = 0
            self.params['teff_2'] = 0
            self.params['logg_1'] = 0
            self.params['logg_2'] = 0
            self.params['logl_1'] = 0
            self.params['logl_2'] = 0
            # self.params['f_contr'] = 0.5

        if self.interpolate_flux:
            f_1 = 10 ** self.params['logl_1']
            f_2 = 10 ** self.params['logl_2']
            flux_ratio = f_1 / (f_1 + f_2)
            self.params['f_contr'] = flux_ratio

            # else:
            #     print("Isochrone labels not found in model labels. Please add 'mass', 'age', and 'fe_h' to model labels for interpolation")

        # elapsed_time = time.time() - start_time
        # print(f"Time taken for interpolation computation: {elapsed_time:.4f} seconds")


        if return_values:
            return interpolate_1, interpolate_2

    # For single parameter retrieval (E.g. teff_1 not teff)
    def get_param(self, param):
        if param[:-2] in ['teff', 'logg', 'logg'] and self.interpolator is not None:
            # We are trying to get a value we need to interpolate for.
            self.interpolate()
            return self.params[param]
        else:
            return self.params[param]

    def get_labels(self):
        return [label for label in self.model_labels.values()]
    
    def get_unique_labels(self):
        return self.unique_labels
    
    def get_component_labels(self, component):
        return [label for label in self.model_labels.values() if label.split('_')[-1] == str(component)]

    def get_params(self, values_only=False, exclude_fixed=False):


        if values_only:
            if exclude_fixed:
                mask = [key.split('_')[0] not in self.fixed_labels for key in self.params.keys()]
                return np.array([float(param) for param in self.params.values()])[mask]
            else:
                return np.array([float(param) for param in self.params.values()])
        else:
            if exclude_fixed:
                return {key: value for key, value in self.params.items() if key.split('_')[0] not in self.fixed_labels}
            else:
                # params_dict = {x: y for x, y in self.params.items()}
                return self.params
    
    # Sets a parameter for all components at once. E.g. set rv bounds for rv_1 and rv_2 simultaneously
    def set_bounds(self, param, bounds=(-np.inf, np.inf)):
        if param in self.unique_labels:
            for i in range(self.components):
                self.bounds[param + '_' + str(i+1)] = bounds
        elif param in self.single_labels:
            self.bounds[param] = bounds
        elif param in self.derived_labels:
            # This is a new parameter that is not fit but is an output.
            # Assume it's for each component and not unique.
            for i in range(self.components):
                print("Setting bound on: ", param + '_' + str(i+1), " to ", bounds)
                self.bounds[param + '_' + str(i+1)] = bounds

    
    # Returns initial parameters as a dictionary without suffixes for each component. E.g. fe_h_1 = 1 -> fe_h = 1
    def get_component_params(self, component, values_only=False, exclude=[]):
        params = []
        for label in self.get_component_labels(component):
            if label[:-2] not in exclude:
                params.append(float(self.params[label]))

        if values_only:
            return params
        else:
            params_dict = {}
            i = 0
            for label in self.get_unique_labels():
                if label not in exclude:
                    params_dict[label] = params[i]
                    i += 1

            return params_dict
    
    # Returns the index of a label in the model_labels dictionary (Enum-ify)
    def label(self, label, comp=None):
        # E.g. if label = 'rv' and comp = 1, return index of rv_1
        if comp is not None:
            return self.indices[label + '_' + str(comp)]
        else:
            return self.indices[label]
    
    def get_comp_mask(self, component):
        # Return a mask for labels in each component
        return [label[-2:] == '_' + str(component) for label in self.model_labels]


    # Returns bounds as an array formatted for curve_fit
    def get_bounds(self, type='list'):
        if type == 'list':
            # Get first bound from each item in dictionary
            bounds_lower = [float(bound[0]) for bound in self.bounds.values()]
            bounds_upper = [float(bound[1]) for bound in self.bounds.values()]

            # Return as tuple instead of list of lists.
            return [tuple(bounds_lower), tuple(bounds_upper)]
        else:
            bounds = [(float(bound[0]), float(bound[1])) for bound in self.bounds.values()]
            return bounds
    
    # Sets a parameter for all components at once. E.g. set rv paramater value for rv_1 and rv_2 simultaneously
    def set_param(self, param, value):
        if param in self.unique_labels:
            for i in range(self.components):
                self.params[param + '_' + str(i+1)] = value
        # New parameter not in model_labels
        elif param in self.model_labels:
            self.params[param] = value
        else:
            print("Parameter " + param +  " not found in model labels. If this is a new parameter or fixed parameter, add it using add_param() first.")

    # Sets all parameters of the model at once
    def set_params(self, params):
        # print("TYPE ", type(params))
        # print("PARAMS: ", params)

        if type(params) is dict:
            for key, value in params.items():
                self.params[key] = value
        else:
            # This is an array of values.
            #  Model labels we are trying to fit
            fit_labels = [label for label in list(self.model_labels) if label.split('_')[0] not in self.fixed_labels]

            if len(fit_labels) == len(params):
                for i, label in enumerate(fit_labels):
                    self.params[label] = params[i]
            else:
                print(fit_labels, len(fit_labels))
                print(fit_labels, params)
                raise ValueError("Error: trying to set all parameters at once but the number of parameters does not match the number of labels in the model.")

    def add_param(self, param, value):
        self.unique_labels.append(param)
        for i in range(self.components):
            self.model_labels[param + '_' + str(i+1)] = param + '_' + str(i+1)
            self.params[param + '_' + str(i+1)] = value

            if param not in self.fixed_labels:
                self.bounds[param + '_' + str(i+1)] = (-1e10, 1e10)

    def generate_model(self):
        self.wavelengths, self.flux, sigma, self.model_flux, unmasked_iter1 = af.return_wave_data_sigma_model(self, self.spectrum, self.same_fe_h) 
        
    def get_residual(self):
        clipped_flux = np.clip(self.flux, 0, 1)
        r = 100 * np.sum(abs(self.model_flux - clipped_flux)) / len(clipped_flux)

        return r
    
    # def get_rchi2(self, with_bounds=True, clipped=True):
    #     # Manual bound checking for interpolated values
    #     if clipped:
    #         clipped_flux = np.clip(self.flux, 0, 1)
    #     else:
    #         clipped_flux = self.flux

    #     r = np.sum((self.model_flux - clipped_flux) ** 2) / (len(clipped_flux) - len(self.params))

    #     # Manual band enforcement for all parameters
    #     if self.unnormalized_bounds and with_bounds:
    #         # print(self.unnormalized_bounds)
    #         for param, bounds in zip(self.get_params(), self.unnormalized_bounds):
    #             # print(f"Checking parameter {param} is outside bounds {self.unnormalized_bounds[param]} with value {self.params[param]}. Penalising residual. {self.unnormalized_bounds[param][1]}")
    #             if self.params[param] > self.unnormalized_bounds[param][1] or self.params[param] < self.unnormalized_bounds[param][0]:
    #                 # print(f"Parameter {param} is outside bounds {self.unnormalized_bounds[param]} with value {self.params[param]}. Penalising residual.")
    #                 r = r * 1e10

    #             # Prevent the RVs from being too close together. This can force bad results be careful.
    #         # if abs(self.params['rv_1'] - self.params['rv_2']) < 60:
    #         #     r = r * 1e10

    #     return r
    
    # New implementation of rchi2 via Perplexity

    def get_rchi2(self, with_bounds=True, clipped=True):
        # Clip flux values if requested
        if clipped == True:
            clipped_flux = np.clip(self.flux, 0, 1.001)
        # Else check if clipped is numerical
        elif isinstance(clipped, (int, float)):
            clipped_flux = np.clip(self.flux, 0, clipped)
        else:
            clipped_flux = self.flux

        # # Calculate residuals (model - observed)
        # residuals = self.model_flux - clipped_flux

        # # Check for valid uncertainties
        # if not hasattr(self, 'flux_uncertainty') or self.flux_uncertainty is None:
        #     raise ValueError("Uncertainty array must be defined for χ² calculation")

        # # Mask points with invalid uncertainties (<=0)
        # valid_mask = self.flux_uncertainty > 0
        # residuals = residuals[valid_mask]
        # uncertainty = np.sqrt(self.flux_uncertainty[valid_mask])

        # # Compute χ² and reduced χ²
        # chi2 = np.sum((residuals / uncertainty) ** 2)
        # n_data = len(residuals)
        # n_params = len(self.params)
        # dof = n_data - n_params
        
        # # Degrees of freedom (DoF) should be positive
        # if dof <= 0:
        #     return np.inf  # Avoid division by zero/negative DoF
        # rchi2 = chi2 / dof


        # rchi2_single = np.median(np.abs(self.flux - self.model_flux_single)/ self.flux_uncertainty)
        # rchi2 = np.median(np.abs(self.flux - self.model_flux)/ self.flux_uncertainty)

        # rchi2 = np.median(np.abs(clipped_flux - self.model_flux)/ np.sqrt(self.flux_uncertainty)) # Take sqrt if using model unc
        # Print length of each array
        # print("Length of clipped flux: ", len(clipped_flux))
        # print("Length of model flux: ", len(self.model_flux))
        # print("Length of flux uncertainty: ", len(self.flux_uncertainty))
        # print("Length of interp ", len(np.interp(self.wavelengths, self.wavelengths_single, self.flux_uncertainty)))
        rchi2 = np.median(np.abs(clipped_flux - self.model_flux)/ self.flux_uncertainty)


        # Manual band enforcement for all parameters
        if self.unnormalized_bounds and with_bounds:
            # print(self.unnormalized_bounds)
            for param, bounds in zip(self.get_params(), self.unnormalized_bounds):
                # print(f"Checking parameter {param} is outside bounds {self.unnormalized_bounds[param]} with value {self.params[param]}. Penalising residual. {self.unnormalized_bounds[param][1]}")
                if self.params[param] > self.unnormalized_bounds[param][1] or self.params[param] < self.unnormalized_bounds[param][0]:
                    # print(f"Parameter {param} is outside bounds {self.unnormalized_bounds[param]} with value {self.params[param]}. Penalising residual.")
                    rchi2 *= 1e10

        return rchi2

    def plot(self, title_text="", lines=None, line_buffer=3, no_lines=5, random_lines=False, component_fluxes=False, component_offset=0, vlines=True, show_plot=True):
        global important_lines
        reg_buf = line_buffer
        # Initialize lists to hold legend handles and labels
        handles, labels = [], []



        if lines is None:
            if random_lines:
                # Choose 10 random indices between 0 and the number of important lines
                random_indices = np.random.choice(len(important_lines), no_lines, replace=False)
                # Get the corresponding lines
                plot_lines = [important_lines[i] for i in random_indices]
            else:
                plot_lines = important_lines[0:no_lines]
                # Remove line at index 2
                plot_lines.pop(2)
                plot_lines.append(important_lines[no_lines+1])

            no_plots = no_lines # How many lines to show from important lines
        else:
            plot_lines = lines
            no_plots = len(lines)

        if self.wavelengths.size > 0 and self.flux.size > 0 and self.model_flux.size > 0:
            # Clip the value to max 30
            figsize = min(60, no_plots * 5)
            fig, axes = plt.subplots(1, no_plots, figsize=(figsize, 5), sharey=True)
            # Iterate over each line and corresponding subplot
            for i, line in enumerate(plot_lines[0:no_plots]):



                # Define the region to plot: line ± 5 Å
                if lines is None:
                    line_wvl = line[0]
                else:
                    line_wvl = line

                    # Check if the line is in the important lines list and use this precise wavelength instead of the one passed in manually.
                    for important_line in important_lines:
                        wavelength, name1, name2 = important_line
                        if abs(line - wavelength) <= 1:
                            line = wavelength
                            line_wvl = wavelength

                if "H" in important_lines[i][1] or "H" in important_lines[i][2]:
                    line_buffer = 10
                else:
                    line_buffer = reg_buf

                min_wave = line_wvl - line_buffer
                max_wave = line_wvl + line_buffer
                
                # Select data within the specified wavelength range
                mask = (self.wavelengths >= min_wave) & (self.wavelengths <= max_wave)
                
                # Plot data and model in the corresponding subplot
                h1, = axes[i].plot(self.wavelengths[mask], self.flux[mask], label='Observed Data')
                h2, = axes[i].plot(self.wavelengths[mask], self.model_flux[mask], label='Model Fit', linestyle='--')

                

                # Model components
                # h4, = axes[i].plot(self.wavelengths[mask], self.model_flux[mask] * self.params['f_contr'], label='Model Fit Primary', linestyle='--')
                # h5, = axes[i].plot(self.wavelengths[mask], self.model_flux[mask] * (1 - self.params['f_contr']), label='Model Fit Primary', linestyle='--')
                if component_fluxes:
                    h4, = axes[i].plot(self.wavelengths[mask], self.component_model_fluxes["model_component_1"][mask], label='Model Fit Primary', linestyle='--', c='r')
                    h5, = axes[i].plot(self.wavelengths[mask], self.component_model_fluxes["model_component_2"][mask], label='Model Fit Primary', linestyle='--', c='b')

                # h6, = axes[i].plot(self.wavelengths[mask], self.component_model_fluxes["model_component_1"][mask] + self.component_model_fluxes["model_component_2"][mask], label='Model Fit Primary', linestyle='--', c='g')

                difference = abs(self.model_flux[mask] - self.flux[mask])
                h3, = axes[i].plot(self.wavelengths[mask], difference, label='Model Delta', linestyle='--')
                axes[i].fill_between(self.wavelengths[mask], 0, difference, color='gray', alpha=0.3)

                clipped_flux = np.clip(self.flux, 0, 1)
                h6, = axes[i].plot(self.wavelengths[mask], clipped_flux[mask], label='Clipped Data', linestyle='--', alpha=0.3)

                # axes[i].set_ylim(np.min(self.flux[mask]) - 0.1, 1.2)

                # vlines = [7771.94, 7774.17, 7775.39]
                # for vline in vlines:
                #     axes[i].vlines(vline, 0, 1, colors='r', linestyles='dotted')
                #     axes[i].vlines(vline + 4, 0, 1, colors='orange', linestyles='dotted')

                # Plot vertical line at the line wavelength
                if vlines:
                    shifted_line = af.rv_shift(self.params['rv_1'] - self.params['rv_1'], line_wvl)
                    axes[i].vlines(shifted_line, 0, 1.15, colors='r', linestyles='dotted')
                    # Annotate the vlines with the radial velocity shifts
                    axes[i].text(shifted_line + 0.2, 1.13, f'$\\lambda_1$: {round(shifted_line, 3)} Å', color='r', bbox=dict(facecolor='white', edgecolor='none', boxstyle='round, pad=0.3'))

                    shifted_line = af.rv_shift(self.params['rv_1'] - self.params['rv_2'], line_wvl)
                    axes[i].vlines(shifted_line, 0, 1.05, colors='blue', linestyles='dotted')
                    # axes[i].text(shifted_line + 0.2, 1.03, f'RV2: {round(self.params["rv_2"], 0)} km/s', color='blue', bbox=dict(facecolor='white', edgecolor='none', boxstyle='round, pad=0.3'))
                    axes[i].text(shifted_line + 0.2, 1.03, f'$\\lambda_2$: {round(shifted_line, 3)} Å', color='blue', bbox=dict(facecolor='white', edgecolor='none', boxstyle='round, pad=0.3'))
                    

                # Set subplot title and labels
                skip_title = False
                if lines is None:
                    axes[i].set_title(f'{line[1]} ({line[0]} Å)')
                else:
                    # Check if any of the plot lines are in the important lines list and use the name instead of the wavelength
                    for important_line in important_lines:
                        wavelength, name1, name2 = important_line
                        if abs(line - wavelength) <= 1:
                            axes[i].set_title(f'{name1} ({line} Å)')
                            skip_title = True

                    if not skip_title:
                        axes[i].set_title(f'({line} Å)')


                axes[i].set_xlabel('Wavelength (Å)')
                x_ticks = np.arange(min_wave, max_wave + 1, line_buffer / 2)
                x_ticks = [int(x) for x in x_ticks]
                axes[i].set_xticks(x_ticks)

                if i == 0:
                    axes[i].set_ylabel('Flux')
                    handles.extend([h1, h2, h3])
                    labels.extend(['Observed Data', 'Binary Model', 'Model Residual'])

                if component_offset < 0:
                    component_offset = 0
                axes[i].set_ylim(-0.1, 1.2 + component_offset)

                if i == 0:
                    # Draw a vertical line on the plot from y = 0.2 to y = 0.8 at x = line_wvl
                    axes[i].axvline(x=line_wvl - line_buffer - 2, ymin=0.2, ymax=0.8, color='black', linestyle='-')
                    # Add a scatter point corresponding to the flux ratio on the line
                    axes[i].scatter(line_wvl - line_buffer - 2, self.params["f_contr"], color='r', s=25, marker='o')
                    axes[i].scatter(line_wvl - line_buffer - 2, self.params["f_contr"], color='r', s=255, marker='_')

            # Adjust layout to prevent overlap
            model_agreement_percentage = 100 * np.sum(abs(self.model_flux - self.flux)) / len(self.flux)
            # residual = np.sum(residuals**2) / (len(residuals) - len(model_parameters))

            # model_agreement_percentage = model_agreement_percentage ** 2
            model_rchi = self.get_rchi2()
            # model_rchi = 1e3 * model_rchi
            
            if self.interpolator is not None:
                title = "ID: " + str(self.id) + \
                    "   Agreement: " + str(model_agreement_percentage) + "%" \
                    "   $\\chi^{2}$: " + str(round(model_rchi, 6))
            else:
                title = "ID: " + str(self.id) + \
                    "   Agreement: " + str(round(model_agreement_percentage, 4)) + "%" \
                    "   $\\chi^{2}$: " + str(round(model_rchi, 6))
            
            title = title + " " + title_text
            plt.suptitle(title)
            # One legend for all plots, positioned in the center underneath
            # axes[0].legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=no_plots)
            fig.legend(handles=handles, labels=labels, loc='upper center', bbox_to_anchor=(0.5, -0.01), ncol=no_plots)


            # For each subplot, set the y-axis to be determined by the min and max of the flux values

            # Annotate the figure
            x_start = 0.4
            if self.interpolate_flux:
                annotation1 = \
                        "   rv$_1$: " + str(round(self.params["rv_1"], 1)) + \
                        "   teff$_1$: " + str(round(self.params["teff_1"] * 1e3)) + \
                        "   logg$_1$: " + str(round(self.params["logg_1"], 3)) + \
                        "   mass$_1$: " + str(round(self.params["mass_1"], 3))
                        # "   age$_1$: " + str(round(self.params["age_1"], 4)) + \
                        # "   metallicity$_1$: " + str(round(self.params["FeH_1"], 4))
                
                annotation2 = \
                        "   rv$_2$: " + str(round(self.params["rv_2"], 1)) + \
                        "   teff$_2$: " + str(round(self.params["teff_2"] * 1e3)) + \
                        "   logg$_2$: " + str(round(self.params["logg_2"], 3)) + \
                        "   mass$_2$: " + str(round(self.params["mass_2"], 3))
                        # "   age$_2$: " + str(round(self.params["age_2"], 4)) + \
                        # "   metallicity$_2$: " + str(round(self.params["FeH_2"], 4))

                fig.text(x_start - 0.12, -0.2, "Flux Ratio: " + str(round(self.params["f_contr"], 4)), ha='center', va='center', fontsize=18)
                fig.text(x_start + 0.12, -0.15, annotation1, ha='center', va='center', fontsize=18)
                fig.text(x_start + 0.12, -0.25, annotation2, ha='center', va='center', fontsize=18)
            fig.text(0.05, 0, f'$RV_1$: {round(self.params["rv_1"], )} km/s', color='r', bbox=dict(facecolor='white', edgecolor='none', boxstyle='round, pad=0.3'))
            fig.text(0.05, -0.05, f'$RV_2$: {round(self.params["rv_2"], )} km/s', color='blue', bbox=dict(facecolor='white', edgecolor='none', boxstyle='round, pad=0.3'))

            plt.tight_layout()

            if show_plot:
                plt.show()  
            else:
                plt.savefig(f'plots/2/spec_{title_text}.png')
                plt.close()

                # params = self.get_params(values_only=True)
                # params_list = ', '.join(map(str, params))
                # with open("animations/fit_results.txt", "a") as f:
                #     f.write(f"{self.id}, {self.get_residual()}, {self.get_rchi2()}, {params_list}\n")


        else:
            print('No data to plot')

    # Plot a specific wavelength region of the spectra including the model, observed flux, and residual
    def plot_window(self, min_wave, max_wave, title_text="", show_plot=True):
        if self.wavelengths.size > 0 and self.flux.size > 0 and self.model_flux.size > 0:
            # Select data within the specified wavelength range
            mask = (self.wavelengths >= min_wave) & (self.wavelengths <= max_wave)

            if not np.any(mask):
                print("No data in the specified wavelength range.")
                return

            # Plot data
            fig, ax = plt.subplots(figsize=(10, 5))
            ax.scatter(self.wavelengths[mask], self.flux[mask], label='Observed Data', color='blue', s=2)
            ax.plot(self.wavelengths[mask], self.flux[mask], lw=2, color='C0', alpha=0.7)
            ax.scatter(self.wavelengths[mask], self.model_flux[mask], label='Model Fit', linestyle='--', color='orange', s=2)
            residual = self.flux[mask] - self.model_flux[mask]
            ax.plot(self.wavelengths[mask], residual, label='Residual', linestyle='--', color='green', lw=2)
            ax.fill_between(self.wavelengths[mask], 0, residual, color='gray', alpha=0.3)

            # Add labels and title
            ax.set_xlabel('Wavelength (Å)')
            ax.set_ylabel('Flux')
            ax.set_title(f'Spectral Window Plot {title_text}')
            ax.legend()

            # Show or save the plot
            if show_plot:
                plt.show()
            else:
                plt.savefig(f'plots/2/window_{title_text}.png')
                plt.close()
        else:
            print('No data to plot')

    def plot_residual(self, title_text="", show_plot=True):
        if self.wavelengths.size > 0 and self.flux.size > 0 and self.model_flux.size > 0:
            fig, ax = plt.subplots(figsize=(20, 5))
            ax.plot(self.wavelengths, self.flux - self.model_flux, label='Residual')
            ax.fill_between(self.wavelengths, 0, self.flux - self.model_flux, color='gray', alpha=0.3)
            ax.set_xlabel('Wavelength (Å)')
            ax.set_ylabel('Residual')
            ax.set_title('Residual Plot ' + title_text)
            plt.legend()
            if show_plot:
                plt.show()
            else:
                plt.savefig(f'plots/2/residual_{title_text}.png')
                plt.close()
        else:
            print('No data to plot')

    def plot_hr(self, background_data):
        GALAH_data = background_data

        plt.figure(figsize=(15, 7))
        cleaned_data = GALAH_data[['teff', 'logg']].replace([np.inf, -np.inf], np.nan).dropna()

        # Plot the 2D histogram
        plt.hexbin(cleaned_data['teff'], cleaned_data['logg'], gridsize=100, cmap='gray_r', alpha=0.5)


        # Select the sample stars in the galah data and show in blue
        unresolved = GALAH_data[GALAH_data['sobject_id'] == self.id]
        plt.scatter(unresolved['teff'], unresolved['logg'], alpha=0.5, color='blue', s=25, label='Unresolved')
        
        # plt.scatter(sample['teff_1'] * 1e3, sample['logg_1'], alpha=0.5, color='red', s=1, label='Primary')
        # plt.scatter(sample['teff_2'] * 1e3, sample['logg_2'], alpha=0.5, color='green', s=1, label='Secondary')

        # plt.

        plt.rcParams.update({'font.size': 20})
        plt.xlabel('teff')
        plt.ylabel('logg')
        plt.ylim(-1, 5)
        plt.legend()
        plt.gca().invert_xaxis()
        plt.gca().invert_yaxis()
        plt.xlim(8000, 3000)
        plt.tight_layout()

        plt.show()

    def plot_binned_residuals(self, bin_size=60):
        bin_size = bin_size  # Bin width in Angstroms

        # Get global wavelength range for this model
        global_bins = np.arange(np.nanmin(self.wavelengths), np.nanmax(self.wavelengths), bin_size)
        bin_centers = 0.5 * (global_bins[:-1] + global_bins[1:])  # Center of each bin

        # Compute residuals
        residuals = self.flux - self.model_flux

        # Mask regions with no data
        valid_mask = ~np.isnan(self.flux) & ~np.isnan(self.model_flux)
        wave_valid = self.wavelengths[valid_mask]
        residuals_valid = residuals[valid_mask]

        # Bin residuals using the global binning scheme
        binned_residuals = np.full(len(bin_centers), np.nan)  # NaN for missing bins

        for i in range(len(global_bins) - 1):
            mask = (wave_valid >= global_bins[i]) & (wave_valid < global_bins[i+1])
            if np.any(mask):  # Only store bins that have data
                binned_residuals[i] = np.nanmean(np.abs(residuals_valid[mask]))  # Mean absolute residual

        # Plot results
        plt.figure(figsize=(10, 5))
        plt.bar(bin_centers, binned_residuals, width=bin_size, alpha=0.5, label=f'Binned |Residuals| ({bin_size} Å)', color='blue')
        plt.plot(bin_centers, binned_residuals, 'o-', label='Scatter Points', color='r')
        plt.axhline(0, color='gray', linestyle='--')
        plt.xlabel("Wavelength (Å)")
        plt.ylabel("Residuals (Flux - Model Flux)")
        plt.legend()
        plt.title("Binned Residuals for Model")

        # Set x-axis limits to valid data range
        plt.xlim(global_bins[0], global_bins[-1])

        plt.show()

        return bin_centers, binned_residuals

    def plot_spectrum(self, mask=None, plot_width=None, show_lines=False):
        # Define the wavelength mask if not provided
        if mask is None:
            mask = (self.wavelengths >= 5708) & (self.wavelengths <= 5718)
        else:
            mask = (self.wavelengths >= mask[0]) & (self.wavelengths <= mask[1])

        masked_wavelengths = self.wavelengths[mask]
        masked_flux = self.flux[mask]
        masked_flux_uncertainty = self.flux_uncertainty[mask]

        if len(masked_wavelengths) == 0:
            print("No data in the specified wavelength range.")
            return

        plt.rcParams.update({'font.size': 20})
        fig_size = min(2 * (max(masked_wavelengths) - min(masked_wavelengths)), 20)
        fig_size = max(20, fig_size)
        if plot_width is not None:
            fig_size = plot_width

        plt.figure(figsize=(fig_size, 4))
        plt.title(f"Observed Spectrum for {self.id}")

        # Plot the spectrum with error bars
        plt.fill_between(masked_wavelengths, masked_flux - masked_flux_uncertainty, masked_flux + masked_flux_uncertainty, color='gray', alpha=0.1)
        plt.errorbar(masked_wavelengths, masked_flux, yerr=masked_flux_uncertainty, fmt='none', ecolor='black', elinewidth=2, capsize=3, alpha=0.4)
        plt.plot(masked_wavelengths, masked_flux, label='Observed Spectrum', color='black', lw=2, alpha=0.8)

        plt.xlabel('Wavelength (Å)')
        plt.ylabel('Flux')
        plt.ylim(min(masked_flux) - 0.2, 1.08)

        if show_lines:
            # Detect lines within 5 angstrom of the center of the mask from important lines
            lines = [line for line in important_lines if abs(line[0] - np.mean(masked_wavelengths)) < 5]

            # For each line, find the closest wavelength in the spectrum and its flux
            lines_with_depth = []
            if len(lines) > 0:
                for line in lines:
                    shifted_line = af.rv_shift(self.params['rv_1'] - self.params['rv_1'], line[0])
                    idx = np.abs(masked_wavelengths - shifted_line).argmin()
                    line_flux = masked_flux[idx]
                    depth = 1 - line_flux  # or continuum - line_flux if not normalized
                    lines_with_depth.append((line, depth))

            if len(lines_with_depth) > 0:
                # Get the line with the largest depth
                deepest_line, max_depth = max(lines_with_depth, key=lambda x: x[1])
                lines = [deepest_line]

                # Plot the lines
                for i, line in enumerate(lines):
                    shifted_line = af.rv_shift(self.params['rv_1'] - self.params['rv_1'], line[0])
                    
                    if min(masked_wavelengths) <= shifted_line <= max(masked_wavelengths):
                        plt.axvline(x=shifted_line, color='darkred', linestyle='--')

                    shifted_line_2 = af.rv_shift(self.params['rv_1'] - self.params['rv_2'], line[0])
                    if min(masked_wavelengths) <= shifted_line_2 <= max(masked_wavelengths):
                        plt.axvline(x=shifted_line_2, color='darkblue', linestyle='--')

                    # Draw a horizontal line between the two lines
                    y_connector = 0.90
                    plt.plot(
                        [shifted_line, shifted_line_2], [y_connector, y_connector],
                        color='black', linestyle='--', lw=1
                    )

                    midpoint_x = (shifted_line + shifted_line_2) / 2
                    # Annotate above the horizontal connector
                    plt.text(
                        midpoint_x, y_connector + 0.05,  # 0.03 above the connector line
                        f'{line[1]} ({round(line[0])} Å)', 
                        color='black', 
                        ha='center',  # Center horizontally
                        va='bottom',  # Align text baseline to bottom of text
                        bbox=dict(facecolor='white', edgecolor='black', boxstyle='round, pad=0.3'),
                        fontsize=10
                    )

                    # Add a vertical line between the midpoint and the textbox
                    plt.plot(
                        [midpoint_x, midpoint_x], [y_connector, y_connector + 0.05],
                        color='black', linestyle='--', lw=1
                    )

        plt.xticks(np.arange(min(masked_wavelengths), max(masked_wavelengths)+1, 4))

        plt.legend(loc='lower right', fontsize=16, frameon=True)
        plt.show()

    # Compare the single star and binary fit plots
    def plot_comparison(self, mask=None, show_lines=False, show_lines_comp=False, plot_width=None):
        if mask is None:
            mask = (self.wavelengths >= 5708) & (self.wavelengths <= 5718)
            # Oxygen triplet
            mask = (self.wavelengths >= 7770) & (self.wavelengths <= 7780)
            # H Alpha
            # mask = (self.wavelengths >= 6552) & (self.wavelengths <= 6572)
        else:
            mask = (self.wavelengths >= mask[0]) & (self.wavelengths <= mask[1])
        
        # same_fe_h = False
        file = '/avatar/buder/GALAH_DR4/analysis_products_single/'+str(self.id)[:6]+'/'+str(self.id)+'/'+str(self.id)+'_single_fit_results.fits'
        global single_results
        try:
            single_results = Table.read(file)
        except:
            print('Single results not available')
            return

        # Apply the mask to the wavelengths and flux
        masked_wavelengths = self.wavelengths[mask]
        masked_flux = self.flux[mask]
        masked_model_flux = self.model_flux[mask]

        # Check if we have any values within the masked wavelengths
        if len(masked_wavelengths) == 0:
            print("No data in the specified wavelength range.")
            return

        # Flux uncertainty may be slightly different in length from the model flux by a few points.
        # Interpolate the flux uncertainty to match the model flux
        masked_flux_uncertainty = self.flux_uncertainty[mask]
        # = np.interp(self.wavelengths, self.wavelengths_single, self.flux_uncertainty)

        # rchi2_single = np.median(np.abs(self.flux - self.model_flux_single)/ self.flux_uncertainty)
        # rchi2 = np.median(np.abs(self.flux - self.model_flux)/ self.flux_uncertainty)

        rchi2_single = self.rchi2_single
        rchi2 = self.get_rchi2()
        


        # Font size increase
        plt.rcParams.update({'font.size': 20})
        fig_size = min(2 * (max(masked_wavelengths) - min(masked_wavelengths)), 20)
        fig_size = max(20, fig_size)

        if plot_width is not None:
            fig_size = plot_width

        plt.figure(figsize=(fig_size, 4))
        plt.title(f"Binary and Single Star Fits for {self.id}")
        # Fill between the error bars
        plt.fill_between(masked_wavelengths, masked_flux - masked_flux_uncertainty, masked_flux + masked_flux_uncertainty, color='gray', alpha=0.1)
        plt.errorbar(masked_wavelengths, masked_flux, yerr=masked_flux_uncertainty, fmt='none', ecolor='black', elinewidth=2, capsize=3, label='Observed Spectra', alpha=0.7)

        plt.plot(masked_wavelengths, masked_model_flux, label=f"Binary {round(self.params['rv_1'])} km$s^{{-1}}$ & {round(self.params['rv_2'])} km$s^{{-1}}$", color='darkblue', lw=2, alpha=0.7)
        plt.plot(masked_wavelengths, self.model_flux_single[mask], label=f"Single {round(single_results['rv_gauss'][0])}  km$s^{{-1}}$", color='darkred', lw=2, alpha=0.7)
        # Annotate the chi2 values in the bottom left corner
        plt.text(0.02, 0.2, f"$\\chi^{2}$: {rchi2:.2f}", transform=plt.gca().transAxes, fontsize=14, verticalalignment='top',color='darkblue')
        plt.text(0.02, 0.1, f"$\\chi^{2}$: {rchi2_single:.2f}", transform=plt.gca().transAxes, fontsize=14, verticalalignment='top', color='darkred')
        plt.xlabel('Wavelength (Å)')
        plt.ylabel('Flux')

        plt.ylim(min(min(masked_flux), min(masked_model_flux), min(self.model_flux_single[mask])) - 0.2, 1.08)
    
        if show_lines:
            # Detect lines within 5 angstrom of the center of the mask from important lines
            lines = [line for line in important_lines if abs(line[0] - np.mean(masked_wavelengths)) < 5]

            # For each line, find the closest wavelength in the spectrum and its flux
            lines_with_depth = []
            if len(lines) > 0:
                for line in lines:
                    shifted_line = af.rv_shift(self.params['rv_1'] - self.params['rv_1'], line[0])
                    idx = np.abs(masked_wavelengths - shifted_line).argmin()
                    line_flux = masked_flux[idx]
                    depth = 1 - line_flux  # or continuum - line_flux if not normalized
                    lines_with_depth.append((line, depth))

            if len(lines_with_depth) > 0:
                # Get the line with the largest depth
                deepest_line, max_depth = max(lines_with_depth, key=lambda x: x[1])
                lines = [deepest_line]

                # Plot the lines
                for i, line in enumerate(lines):
                    shifted_line = af.rv_shift(self.params['rv_1'] - self.params['rv_1'], line[0])
                    
                    if min(masked_wavelengths) <= shifted_line <= max(masked_wavelengths):
                        plt.axvline(x=shifted_line, color='darkred', linestyle='--')
                        # Annotate at 80% of y-axis range
                        # y_pos = 1.05 * (plt.ylim()[1] - plt.ylim()[0]) + plt.ylim()[0]
                        # plt.text(shifted_line + 0.2, 1.05,  s=f'{line[1]} ({round(line[0])} Å)', color='darkred', 
                        #         bbox=dict(facecolor='white', edgecolor='none', boxstyle='round, pad=0.3'))

                    if show_lines_comp == 1:
                        continue
                    shifted_line_2 = af.rv_shift(self.params['rv_1'] - self.params['rv_2'], line[0])
                    if min(masked_wavelengths) <= shifted_line_2 <= max(masked_wavelengths):
                        plt.axvline(x=shifted_line_2, color='darkblue', linestyle='--')
                        # Annotate at 80% of y-axis range
                        
                        # y_pos = 1.01 * (plt.ylim()[1] - plt.ylim()[0]) + plt.ylim()[0]
                        # plt.text(shifted_line_2 + 0.2, 1.03, s=f'{line[2]} ({round(line[0])} Å)', color='darkblue', 
                        #         bbox=dict(facecolor='white', edgecolor='none', boxstyle='round, pad=0.3'))
                        
                    # Draw a horizontal line between the two lines
                    y_connector = plt.ylim()[1] + 0.1 + (i * 0.05) * (plt.ylim()[1] - plt.ylim()[0])
                    plt.plot(
                        [shifted_line, shifted_line_2], [y_connector, y_connector],
                        color='black', linestyle='--', lw=1
                    )

                    midpoint_x = (shifted_line + shifted_line_2) / 2
                    # Annotate above the horizontal connector
                    plt.text(
                        midpoint_x, y_connector + 0.05,  # 0.03 above the connector line
                        f'{line[1]} ({round(line[0])} Å)', 
                        color='black', 
                        ha='center',  # Center horizontally
                        va='bottom',  # Align text baseline to bottom of text
                        bbox=dict(facecolor='white', edgecolor='black', boxstyle='round, pad=0.3'),
                        fontsize=10
                    )

                    # Add a vertical line between the midpoint and the textbox
                    plt.plot(
                        [midpoint_x, midpoint_x], [y_connector, y_connector + 0.05],
                        color='black', linestyle='--', lw=1
                    )

            plt.ylim(plt.ylim()[0], plt.ylim()[1] + 0.5)

        # Annotate the flux ratio
        plt.text(0.02, 0.3, f"Flux Ratio: {round(self.params['f_contr'], 2)}", transform=plt.gca().transAxes, fontsize=16, verticalalignment='top', color='black', alpha=0.9)
        plt.legend(loc='lower right', fontsize=16, frameon=True)
        # plt.savefig(f'comparison_{self.id}.pdf', bbox_inches='tight', dpi=300)
        plt.show()


    def CrossCorrelate(self, shift_range=np.arange(-1000,1000),  normalise=False):

        rv_1 = self.params['rv_1']

        # Define the range of shifts
        shift_range = np.arange(-1200, 1200).astype(np.float64)

        # wavelengths = np.array(self.wavelengths_single) #+ rv_1 # Wavelengths
        # flux_model = self.model_flux_single
        # flux_obs = self.flux

        # # Initialize an array to store cross-correlation results for each shift
        cross_corr_results = np.zeros_like(shift_range, dtype=float)

        deltas = np.zeros_like(shift_range, dtype=float)
        # # Step 3: Calculate the cross-correlation for each shift
        # for i, shift in enumerate(shift_range):
        #     # Translate the model spectrum by the current shift
        #     shifted_wavelengths = af.rv_shift(shift, wavelengths)

        #     # # Calculate the flux at the new RV value
        #     # flux_at_shifted = np.interp(wavelengths, shifted_wavelengths, self.model_flux_single_original)

        #     # # Ensure we only consider where flux_obs and flux_at_shifted overlap. flux_obs can be shorter than flux_at_shifted
        #     # mask = np.isfinite(flux_obs) & np.isfinite(flux_at_shifted)
        #     # flux_at_shifted = flux_at_shifted[0:len(flux_obs)]
            
        #     # # Calculate sum difference between observed and shifted modelled flux
        #     # delta = 1/np.sum((flux_obs - flux_at_shifted) ** 2)
        #     # deltas[i] = delta


        wavelengths_obs = self.wavelengths  # grid for observed flux
        flux_obs = self.flux

        wavelengths_model = self.wavelengths_single  # grid for the model
        model_flux = self.model_flux_single_original

        for i, shift in enumerate(shift_range):
            shifted_wavelengths = af.rv_shift(shift, wavelengths_model)
            # Interpolate model onto observed grid
            flux_at_shifted = np.interp(wavelengths_obs, shifted_wavelengths, model_flux)
            # Now both arrays are on the observed grid
            mask = np.isfinite(flux_obs) & np.isfinite(flux_at_shifted)
            flux_obs_masked = flux_obs[mask]
            flux_at_shifted_masked = flux_at_shifted[mask]


            # Calculate sum difference between observed and shifted modelled flux
            delta = 1/np.sum((flux_obs - flux_at_shifted) ** 2)
            deltas[i] = delta


        # Step 4: Find the best fit shift
        if normalise:
            deltas = deltas / np.max(deltas)
            
        # TODO fit a gaussian instead for the accurate RVs.jpg
        best_fit_index = np.argmax(deltas)
        highest_peak = deltas[best_fit_index]
        # best_fit_shift = shift_range[best_fit_index]


        # print(deltas)
        # find_peaks returns the indices of the peaks.
        peaks, _ = find_peaks(deltas, prominence=highest_peak * 0.1, height=highest_peak * 0.1)  # Negate deltas as find_peaks looks for valleys
        peaks_d1, _ = find_peaks(np.gradient(np.abs(deltas)), prominence=highest_peak * 0.1, height=highest_peak * 0.1)  # Negate deltas as find_peaks looks for valleys
        peaks_d2, _ = find_peaks(np.gradient(np.gradient(np.abs(deltas))), prominence=highest_peak * 0.1, height=highest_peak * 0.1)  # Negate deltas as find_peaks looks for valleys
        peaks_d3, _ = find_peaks(np.gradient(np.gradient(np.gradient(np.abs(deltas)))), prominence=highest_peak * 0.1, height=highest_peak * 0.1)  # Negate deltas as find_peaks looks for valleys

        global single_results
        # The model has already been shifted by the original RV found for this star to match the data. (A shift of 0 will be the best fit, as the model data has already been shifted when imported here).
        # We need to shift the model by the original RV found, to determine the absolute RV for two stars, otherwise we get relative shifts
        shift_range += single_results['rv_gauss'][0] 


        peak_array = [peaks, peaks_d1, peaks_d2, peaks_d3]

        RVs = shift_range
        # After finding peaks:
        absolute_RVs = rv_1 + shift_range[peaks] - shift_range[peaks[0]]
        RVs = rv_1 + shift_range - shift_range[peaks[0]]


        annotate_fs = 11
        plt.figure(figsize=(20, 4))
        plt.plot(RVs, deltas, label='Cross-Correlation Function', c='black')  # Optionally plot in absolute RV
        plt.scatter(absolute_RVs, deltas[peaks], marker="x", color='red', s=150, label='CCF Peaks')
        for i, peak in enumerate(peaks):
            plt.axvline(x=absolute_RVs[i], c='red', ls='dashed', label=f'CCF RV$_{i+1}$', zorder=1)
            plt.annotate(f'CCF RV$_{i+1}$: {absolute_RVs[i]:.0f}', (absolute_RVs[i] + 2, np.max(deltas) * 1.3 * 0.75), textcoords="offset points", xytext=(0,10), ha='center', fontsize=annotate_fs,
            bbox=dict(facecolor='white', edgecolor='none', boxstyle='round,pad=0.5', alpha=0.6), rotation=90)

        plt.axvline(x=single_results['rv_gauss'][0], c='black', ls='dotted', label='GALAH RV$_1$')
        plt.annotate(f'GALAH RV$_1$: {single_results["rv_gauss"][0]:.0f}', (single_results['rv_gauss'][0] + 2, 0), textcoords="offset points", xytext=(0,10), ha='center', fontsize=annotate_fs,
            bbox=dict(facecolor='white', edgecolor='none', boxstyle='round,pad=0.5', alpha=0.6), rotation=90)

        for i, PSO_rv in enumerate([self.params['rv_1'], self.params['rv_2']]):
            if abs(PSO_rv - absolute_RVs[0]) > 10:
                plt.axvline(x=PSO_rv, c='darkred', ls='dashdot', label=f'PSO RV$_{i+1}$')
            else:
                plt.axvline(x=PSO_rv, c='darkred', ls='dashdot', label=f'PSO RV$_{i+1}$', zorder=0)

            plt.annotate(f'RV$_{i+1}$: {PSO_rv:.0f}', (PSO_rv + 2, 0), textcoords="offset points", xytext=(0,10), ha='center', fontsize=annotate_fs,
            bbox=dict(facecolor='white', edgecolor='none', boxstyle='round,pad=0.5', alpha=0.7), rotation=90)

        plt.xlabel('Radial Velocity (km/s)')
        plt.ylabel('Cross-Correlation Value')
        plt.title('Cross-Correlation Function and Radial Velocity Peaks')
        plt.legend(loc='upper right', fontsize=12)

        plt.xlim(-250 + rv_1, 250 + rv_1)

        min_peak = np.min([self.params['rv_1'], self.params['rv_2']])
        max_peak = np.max([self.params['rv_1'], self.params['rv_2']])
        plt.xlim(min_peak - 50, max_peak + 50)

        plt.ylim(0, np.max(deltas) * 1.6)
        plt.savefig(f'cross_correlation_{self.id}.pdf', bbox_inches='tight', dpi=300)
        plt.show()

        return peaks, shift_range, deltas