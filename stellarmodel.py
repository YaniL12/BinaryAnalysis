import numpy as np
import matplotlib.pyplot as plt
import AnalysisFunctions as af
important_lines, important_molecules = af.load_dr3_lines()

class StellarModel:

    # Constructor with default labels and components. Override for custom labels and components
    def __init__(self, labels = ['rv', 'teff', 'logg', 'fe_h', 'vmic', 'vsini'], components = 2, same_fe_h=False):

    # Dictionaries for labels, bounds, and initial values
    # These should be instance level variables not class level. Otherwise they will be shared between all instances.
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


    def get_labels(self):
        return [label for label in self.model_labels.values()]
    
    def get_unique_labels(self):
        return self.unique_labels
    
    def get_component_labels(self, component):
        return [label for label in self.model_labels.values() if label.split('_')[-1] == str(component)]

    def get_params(self):
        return np.array([float(param) for param in self.params.values()])
    
    # Sets a parameter for all components at once. E.g. set rv bounds for rv_1 and rv_2 simultaneously
    def set_bounds(self, param, bounds=(-np.inf, np.inf)):
        if param in self.unique_labels:
            for i in range(self.components):
                self.bounds[param + '_' + str(i+1)] = bounds

    
    # Returns initial parameters as a dictionary without suffixes for each component. E.g. fe_h_1 = 1 -> fe_h = 1
    def get_component_params(self, component, values_only=False):
        params = []
        for label in self.get_component_labels(component):
            params.append(float(self.params[label]))

        if values_only:
            return params
        else:
            params_dict = {}
            for i, label in enumerate(self.get_unique_labels()):
                params_dict[label] = params[i]

            return params_dict
    
    # Returns the index of a label in the model_labels dictionary (Enumify)
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
    def get_bounds(self):
        # Get first bound from each item in dictionary
        bounds_lower = [float(bound[0]) for bound in self.bounds.values()]
        bounds_upper = [float(bound[1]) for bound in self.bounds.values()]

        # Return as tuple instead of list of lists.
        return [tuple(bounds_lower), tuple(bounds_upper)]
    
    # Sets a parameter for all components at once. E.g. set rv paramater value for rv_1 and rv_2 simultaneously
    def set_param(self, param, value):
        if param in self.unique_labels:
            for i in range(self.components):
                self.params[param + '_' + str(i+1)] = value

    # Sets all parameters of the model at once
    # Would be better to pass a dictionary of parameters instead of a list TODO
    def set_params(self, params):
        for i, label in enumerate(self.model_labels):
            self.params[label] = params[i]

    def generate_model(self, spectrum):
        self.wavelengths, self.flux, sigma2_iter1, self.model_flux, unmasked_iter1 = af.return_wave_data_sigma_model(self, spectrum, same_fe_h = False) 
        

    def plot(self):
        global important_lines
        
        if self.wavelengths.size > 0 and self.flux.size > 0 and self.model_flux.size > 0:
            fig, axes = plt.subplots(1, 10, figsize=(30, 5), sharey=True)
            # Iterate over each line and corresponding subplot
            for i, line in enumerate(important_lines[0:10]):
                # Define the region to plot: line ± 5 Å
                line_wvl = line[0]
                min_wave = line_wvl - 5
                max_wave = line_wvl + 5
                
                # Select data within the specified wavelength range
                mask = (self.wavelengths >= min_wave) & (self.wavelengths <= max_wave)
                
                # Plot data and model in the corresponding subplot
                axes[i].plot(self.wavelengths[mask], self.flux[mask], label='Observed Data')
                axes[i].plot(self.wavelengths[mask], self.model_flux[mask], label='Model Fit', linestyle='--')

                difference = abs(self.model_flux[mask] - self.flux[mask])
                axes[i].plot(self.wavelengths[mask], difference, label='Model Delta', linestyle='--')
                axes[i].fill_between(self.wavelengths[mask], 0, difference, color='gray', alpha=0.3)
                
                # Set subplot title and labels
                axes[i].set_title(f'{line[1]} ({line[0]} Å)')
                axes[i].set_xlabel('Wavelength')
                if i == 0:
                    axes[i].set_ylabel('Flux')
                

            # Adjust layout to prevent overlap
            model_agreement_percentage = 100 * np.sum(abs(self.model_flux - self.flux)) / len(self.flux)
            plt.suptitle(model_agreement_percentage)
            plt.tight_layout()
            plt.show()
        else:
            print('No data to plot')