import numpy as np
import matplotlib.pyplot as plt
import AnalysisFunctions as af
from pandas import DataFrame
important_lines, important_molecules = af.load_dr3_lines()

class StellarModel:
    interpolator = None

    # Constructor with default labels and components. Override for custom labels and components
    def __init__(self, id="No ID", labels = ['rv', 'teff', 'logg', 'fe_h', 'vmic', 'vsini'], fixed_labels=[], components = 2, same_fe_h=False, interpolator=None, interpolate_flux=False):

    # Dictionaries for labels, bounds, and initial values
    # These should be instance level variables not class level. Otherwise they will be shared between all instances.
        self.id = id
        self.components = components
        self.model_labels = {}
        self.unique_labels = []
        self.bounds = {}
        self.params = {}
        self.indices = {}
        self.unique_indices = {}
        self.wavelengths = []
        self.flux = []
        self.model_flux = []

        self.param_data = {}

        # Should we use a simple flux model or determine flux ratios from interpolated luminosity?
        self.interpolate_flux = interpolate_flux

        self.fixed_labels = fixed_labels

        # self.unique_labels.append('f_contr')
        self.unique_labels.extend(labels)

        # Only one instance of f_contr. This will need to be changeed for multiple systems where n_comps > 2
        self.model_labels['f_contr'] = 'f_contr'

        for i, label in enumerate(labels):
            self.unique_indices[label] = i

        for j in range(self.components):
            for i, label in enumerate(labels):
                self.model_labels[label + '_' + str(j+1)] = label + '_' + str(j+1)


        # Instantiate bounds and params dictionaries with defaults
        for i, label in enumerate(self.model_labels):
            self.bounds[label] = (-np.inf, np.inf)
            self.params[label] = 0
            self.indices[label] = i

        if interpolator is not None:
            self.interpolator = interpolator

            if all(label in self.unique_labels for label in ['mass']) and all(label in self.fixed_labels for label in ['age', 'metallicity']):
                self.add_param('teff', 0)
                self.add_param('logg', 0)
                self.add_param('logl', 0)
        
        self.param_data = {key: [] for key in self.params.keys()}
        self.param_data['residual'] = []

    def save_data(self):
        for i, param in enumerate(self.params):
            self.param_data[param].append(self.params[param])

        r = 100 * np.sum(abs(np.array(self.model_flux) - np.array(self.flux))) / len(self.flux)
        self.param_data['residual'].append(r)

    def load_data(self, data):
        # Check if data is a DataFrame
        if type(data) is DataFrame:
            row = data[data['sobject_id'] == self.id]
            for col in data.columns:
                # If col exists in the model.params, set the value
                if col in self.get_labels():
                    self.params[col] = row[col].iloc[0]

            spectrum = af.read_spectrum(self.id)
            self.generate_model(spectrum)

    def interpolate(self):
        if self.interpolator is not None:
            # Accepts mass, log(age), metallicity.
            # Outputs Teff, logg, and log(L) bolometric (flux)

            if all(label in self.unique_labels for label in ['mass']) and all(label in self.fixed_labels for label in ['age', 'metallicity']):
                # Both components will have the same starting values.
                # Provide the log of the age in Gyr for interpolation
                # TODO Change this to be dynamic based on the number of components

                # Star 1
                interpolate_1 = self.interpolator(self.params['mass_1'], np.log10(self.params['age_1'] * 1e9), self.params['metallicity_1'])
                # Star 2
                interpolate_2 = self.interpolator(self.params['mass_2'], np.log10(self.params['age_2'] * 1e9), self.params['metallicity_2'])

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
                if np.isnan(self.params['teff_1']) or np.isnan(self.params['teff_2']):
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

            else:
                print("Isochrone labels not found in model labels. Please add 'mass', 'age', and 'metallicity' to model labels for interpolation")

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
            fit_labels = [label for label in self.model_labels if label.split('_')[0] not in self.fixed_labels]
            
            # print(fit_labels, len(fit_labels))
            # print(fit_labels, self.model_labels, self.fixed_labels)

            if len(fit_labels) == len(params):
                for i, label in enumerate(fit_labels):
                    self.params[label] = params[i]
            else:
                raise ValueError("Error: trying to set all parameters at once but the number of parameters does not match the number of labels in the model.")

    def add_param(self, param, value):
        self.unique_labels.append(param)
        for i in range(self.components):
            self.model_labels[param + '_' + str(i+1)] = param + '_' + str(i+1)
            self.params[param + '_' + str(i+1)] = value

            if param not in self.fixed_labels:
                self.bounds[param + '_' + str(i+1)] = (-1e10, 1e10)

    def generate_model(self, spectrum):
        self.wavelengths, self.flux, sigma2_iter1, self.model_flux, unmasked_iter1 = af.return_wave_data_sigma_model(self, spectrum, same_fe_h = False) 
        
    def get_residual(self):
        return 100 * np.sum(abs(self.model_flux - self.flux)) / len(self.flux)
    
    def get_rchi2(self):
        return np.sum((self.model_flux - self.flux) ** 2) / (len(self.flux) - len(self.params))

    def plot(self, title_text=""):
        global important_lines
        
        # Initialize lists to hold legend handles and labels
        handles, labels = [], []

        no_plots = 10
        if self.wavelengths.size > 0 and self.flux.size > 0 and self.model_flux.size > 0:
            fig, axes = plt.subplots(1, no_plots, figsize=(30, 5), sharey=True)
            # Iterate over each line and corresponding subplot
            for i, line in enumerate(important_lines[0:no_plots]):
                # Define the region to plot: line ± 5 Å
                line_wvl = line[0]
                min_wave = line_wvl - 5
                max_wave = line_wvl + 5
                
                # Select data within the specified wavelength range
                mask = (self.wavelengths >= min_wave) & (self.wavelengths <= max_wave)
                
                # Plot data and model in the corresponding subplot
                h1, = axes[i].plot(self.wavelengths[mask], self.flux[mask], label='Observed Data')
                h2, = axes[i].plot(self.wavelengths[mask], self.model_flux[mask], label='Model Fit', linestyle='--')

                difference = abs(self.model_flux[mask] - self.flux[mask])
                h3, = axes[i].plot(self.wavelengths[mask], difference, label='Model Delta', linestyle='--')
                axes[i].fill_between(self.wavelengths[mask], 0, difference, color='gray', alpha=0.3)
                
                # Set subplot title and labels
                axes[i].set_title(f'{line[1]} ({line[0]} Å)')
                axes[i].set_xlabel('Wavelength')
                if i == 0:
                    axes[i].set_ylabel('Flux')
                    handles.extend([h1, h2, h3])
                    labels.extend(['Observed Data', 'Model Fit', 'Model Delta'])

                axes[i].set_ylim(-0.1, 1.2)

            # Adjust layout to prevent overlap
            model_agreement_percentage = 100 * np.sum(abs(self.model_flux - self.flux)) / len(self.flux)
            # residual = np.sum(residuals**2) / (len(residuals) - len(model_parameters))

            # model_agreement_percentage = model_agreement_percentage ** 2
            model_rchi = self.get_rchi2()
            model_rchi = 1e3 * model_rchi
            
            if self.interpolator is not None:
                title = "ID: " + str(self.id) + \
                    "   Agreement: " + str(model_agreement_percentage) + "%" \
                    "   rchi2" + str(model_rchi)
            else:
                title = "ID: " + str(self.id) + \
                    "   Agreement: " + str(round(model_agreement_percentage, 4)) + "%" \
                    "   $r\\chi^{2} 10^{3}$: " + str(round(model_rchi, 6))
            
            title = title + " " + title_text
            plt.suptitle(title)
            # One legend for all plots, positioned in the center underneath
            # axes[0].legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=no_plots)
            fig.legend(handles=handles, labels=labels, loc='upper center', bbox_to_anchor=(0.5, -0.01), ncol=no_plots)

            # Annotate the figure
            x_start = 0.4
            annotation1 = \
                    "   rv$_1$: " + str(round(self.params["rv_1"], 1)) + \
                    "   teff$_1$: " + str(round(self.params["teff_1"] * 1e3)) + \
                    "   logg$_1$: " + str(round(self.params["logg_1"], 3)) + \
                    "   mass$_1$: " + str(round(self.params["mass_1"], 3)) + \
                    "   age$_1$: " + str(round(self.params["age_1"], 4)) + \
                    "   metallicity$_1$: " + str(round(self.params["metallicity_1"], 4))
            
            annotation2 = \
                    "   rv$_2$: " + str(round(self.params["rv_2"], 1)) + \
                    "   teff$_2$: " + str(round(self.params["teff_2"] * 1e3)) + \
                    "   logg$_2$: " + str(round(self.params["logg_2"], 3)) + \
                    "   mass$_2$: " + str(round(self.params["mass_2"], 3)) + \
                    "   age$_2$: " + str(round(self.params["age_2"], 4)) + \
                    "   metallicity$_2$: " + str(round(self.params["metallicity_2"], 4))

            fig.text(x_start - 0.12, -0.2, "Flux Ratio: " + str(round(self.params["f_contr"], 4)), ha='center', va='center', fontsize=18)
            fig.text(x_start + 0.12, -0.15, annotation1, ha='center', va='center', fontsize=18)
            fig.text(x_start + 0.12, -0.25, annotation2, ha='center', va='center', fontsize=18)

            plt.tight_layout()
            plt.show()
        else:
            print('No data to plot')