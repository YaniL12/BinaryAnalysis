import numpy as np
import matplotlib.pyplot as plt
import AnalysisFunctions as af
from pandas import DataFrame
from pandas import Series as Series
important_lines, important_molecules = af.load_dr3_lines()

class StellarModel:
    interpolator = None

    # Constructor with default labels and components. Override for custom labels and components
    def __init__(self, id="No ID", labels = ['rv', 'teff', 'logg', 'FeH', 'vmic', 'vsini'], fixed_labels=[], single_labels=[], components = 2, same_fe_h=False, interpolator=None, interpolate_flux=False, original_bounds=None):

    # Dictionaries for labels, bounds, and initial values
    # These should be instance level variables not class level. Otherwise they will be shared between all instances.
        self.id = id
        self.components = components
        self.same_fe_h = same_fe_h

        self.model_labels = {}
        self.unique_labels = []
        self.single_labels = single_labels

        self.bounds = {}
        self.original_bounds = original_bounds
        self.params = {}
        self.indices = {}
        self.unique_indices = {}
        self.wavelengths = []
        self.flux = []
        self.model_flux = []
        self.component_model_fluxes = {}

        self.param_data = {}

        self.isochrone_table = None

        # Should we use a simple flux model or determine flux ratios from interpolated luminosity?
        self.interpolate_flux = interpolate_flux
        if interpolate_flux and interpolator is None:
            raise ValueError("You have set interpolate_flux to True but have not provided an interpolator. Please provide an interpolator to use this feature.")

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

        # Add single labels to the model labels. These are those that are the same for both components.
        for i, label in enumerate(single_labels):
            self.model_labels[label] = label

        # Instantiate bounds and params dictionaries with defaults
        for i, label in enumerate(self.model_labels):
            self.bounds[label] = (-np.inf, np.inf)
            self.params[label] = 0
            self.indices[label] = i

        if interpolator is not None:
            self.interpolator = interpolator

            # Interpolator requires mass, age, and metallicity to be present in the model labels
            # Mass will be different for each component, the remaining parameters will be the same for both components (age, FeH).
            if same_fe_h and all(label in self.single_labels for label in ['age', 'FeH']) and all(label in self.unique_labels for label in ['mass']):
                self.add_param('teff', 0)
                self.add_param('logg', 0)
                self.add_param('logl', 0)

            # elif all(label in self.unique_labels for label in ['mass', 'age']) and all(label in self.fixed_labels for label in ['FeH']):
            elif all(label in self.unique_labels for label in ['mass', 'age', 'FeH']):
                self.add_param('teff', 0)
                self.add_param('logg', 0)
                self.add_param('logl', 0)

                # self.interpolate()
        
        self.param_data = {key: [] for key in self.params.keys()}
        self.param_data['residual'] = []

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
                # If col exists in the model.params, set the value
                if col in self.get_labels():
                    self.params[col] = row[col].iloc[0]

            
            spectrum = af.read_spectrum(self.id)
            af.load_neural_network(spectrum)
            self.generate_model(spectrum)

        elif isinstance(data, Series):
            for col in data.keys():
                # If col exists in the model.params, set the value
                if col in self.get_labels():
                    # print(col, data[col])
                    self.params[col] = data[col]
            # print(self.params)
            spectrum = af.read_spectrum(self.id)
            af.load_neural_network(spectrum)
            self.generate_model(spectrum)

    # Call this if we have modified one or more parameters manually. This will update the model flux.
    # Effectively a load data call, without loadiing paramaters.
    def regenerate_spectrum(self):
        spectrum = af.read_spectrum(self.id)
        self.generate_model(spectrum)

    def interpolate(self):
        # Accepts mass, log(age), metallicity.
        # Outputs Teff, logg, and log(L) bolometric (flux)

        if self.interpolator == 'trilinear':

            age_query = np.log10(self.params['age'] * 1e9)
            m_h_query = self.params['FeH']
            mass_query_1 = self.params['mass_1'] - 0.1
            mass_query_2 = self.params['mass_2'] - 0.1

            interpolate_1 = af.interpolate_isochrone(mass_query_1, age_query, m_h_query)
            interpolate_2 = af.interpolate_isochrone(mass_query_2, age_query, m_h_query)


            self.params['teff_1'] = (interpolate_1['teff']) / 1000
            self.params['logg_1'] = interpolate_1['logg']
            self.params['logl_1'] = interpolate_1['logl']

            self.params['teff_2'] = (interpolate_2['teff'])  / 1000
            self.params['logg_2'] = interpolate_2['logg']
            self.params['logl_2'] = interpolate_2['logl']

        elif self.interpolator is not None:

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

    def generate_model(self, spectrum):
        self.wavelengths, self.flux, sigma2_iter1, self.model_flux, unmasked_iter1 = af.return_wave_data_sigma_model(self, spectrum, self.same_fe_h) 
        
    def get_residual(self):
        return 100 * np.sum(abs(self.model_flux - self.flux)) / len(self.flux)
    
    def get_rchi2(self):
        # Manual bound checking for interpolated values
        # print(self.params['f_contr'], self.bounds['f_contr'])
        if self.params['f_contr'] > self.bounds['f_contr'][1] or self.params['f_contr'] < self.bounds['f_contr'][0]:
            return 1e10

        return np.sum((self.model_flux - self.flux) ** 2) / (len(self.flux) - len(self.params))

    def plot(self, title_text="", lines=None, line_buffer=10, no_lines=5, random_lines=False, component_fluxes=False, component_offset=0, vlines=True, show_plot=True):
        global important_lines
        
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
            model_rchi = 1e3 * model_rchi
            
            if self.interpolator is not None:
                title = "ID: " + str(self.id) + \
                    "   Agreement: " + str(model_agreement_percentage) + "%" \
                    "   $r\\chi^{2} 10^{3}$: " + str(round(model_rchi, 6))
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

    