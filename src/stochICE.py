

import os
import shutil

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle

#import stochRIVICE
import stochHECRAS
import random


from scipy.stats import genextreme, gumbel_r, norm, weibull_min




class StochICE_HECRAS:
    """
    Class for managing stochastic simulations in ice-flood modeling 
    using both HECRAS and RIVICE models.
    """
    
    def __init__(self, prjDir: str, batch_ID: str, ras_file: str, geo_file: str, 
                 flow_file: str, wse_file: str, NSims: int, ice_jam_reach: int, stochvars: dict):

        self.prjDir = prjDir
        self.ID = batch_ID
        self.ras_file = os.path.join(self.prjDir, ras_file)
        self.geo_file = os.path.join(self.prjDir, geo_file)
        self.flow_file = os.path.join(self.prjDir, flow_file)
        self.geo_file_temp = os.path.join(self.prjDir, geo_file + '.temp')
        self.wse_map_path = os.path.join(self.prjDir, "MonteCarlo", wse_file)

        # Simulation parameters
        self.NSims = NSims
        self.ice_jam_reach=ice_jam_reach
        self.stochVars = stochvars

        # Initialize variables
        self.xs_data = {}
        self.ice_variables = ['Ice Thickness', 'Ice Mann', 'Ice Specific Gravity', 
                              'Ice Is Channel', 'Ice Is OB', 'Ice Friction Angle', 
                              'Ice Porosity', 'Ice K1', 'Ice Max Mean', 'Ice Cohesion', 
                              'Ice Fixed Mann']

        # Initialize methods
        self.print_header()
        self.setup_monte_carlo_dir()
        self.copy_geofile()
        self.get_XS_ice_params()
        self.get_XS_manning()
        self.get_XS_bank_stations()
        self.get_XS_geometry()
        self.get_XS_main_chan_geometry()
        self.get_bridge_data()


        # Call to check for stochastic variables
        self.check_stochVars()

        # Call to generate stochastic values based on stochvars
        self.generate_stochastic_values()

        self.clr=True
        self.stochHECRAS = stochHECRAS.StochHECRAS(self)
        self.stochHECRAS.preprocess_sims()
        self.stochHECRAS.new_launch_sims()
        self.stochHECRAS.make_ensemble_flood_map()

    def print_header(self):
        """
        Prints the header information for the current batch simulation.

        Args:
            None

        Returns:
            None
        """
        print("\n")
        print("----------------------------------------------------------")
        print("      Running stochICE flood modelling system v0.0.1      ")
        print("----------------------------------------------------------\n")
        print(f"Running HEC-RAS project file: {self.ras_file}")
        print(f"HEC-RAS geometry file: {self.geo_file}")
        print(f"HEC-RAS flow file: {self.flow_file}")
        print(f"Batch ID: {self.ID}")
        print(f"Total simulations: {self.NSims}")

    def setup_monte_carlo_dir(self):
        """
        Sets up directories for Monte Carlo simulations. 
        Creates directories for simulation results, geo files, flow files, and tifs.

        Args:
            None

        Returns:
            None
        """
        self.MC_path = os.path.join(self.prjDir, "MonteCarlo")
        self.tif_path = os.path.join(self.MC_path, "SimulationTifs")
        self.geo_files_path = os.path.join(self.MC_path, "GeoFiles")
        self.flow_files_path = os.path.join(self.MC_path, "FlowFiles")

        os.makedirs(self.MC_path, exist_ok=True)
        os.makedirs(self.tif_path, exist_ok=True)
        os.makedirs(self.geo_files_path, exist_ok=True)
        os.makedirs(self.flow_files_path, exist_ok=True)
        print("Monte Carlo directories set up successfully.")

    def copy_geofile(self):
        """
        Creates a backup of the geometry file before starting the simulation.

        Args:
            None

        Returns:
            None
        """
        file_name, file_ext = os.path.splitext(os.path.basename(self.geo_file))
        geo_file_copy = os.path.join(self.prjDir, f"{file_name}_backup{file_ext}")

        if os.path.exists(geo_file_copy):
            shutil.copyfile(geo_file_copy, self.geo_file)
            print(".g0* file copied from original.")
        else:
            shutil.copyfile(self.geo_file, geo_file_copy)
            print("Copy of original .g0* file made.")


    def get_XS_ice_params(self):
        """
        Extracts ice parameters from cross-sections in the geometry file
        and assigns each cross-section to its respective reach based on 
        the line 'River Reach=' (e.g., 'River Reach=River 1').
    
        The reach number is dynamically extracted from the line.
    
        Args:
            None
    
        Returns:
            None
        """
        self.xs_data = {}
        current_reach = None  # Track the current reach dynamically
        
        with open(self.geo_file, 'r') as f:
            flag = False
            for count, line in enumerate(f):
                # Check if we are entering a new reach (match any line containing "River Reach=")
                if "River Reach=" in line:
                    # Extract the reach number from the line (e.g., 'River Reach=River 1')
                    current_reach = int(line.split("River Reach=River")[1].split()[0].strip(","))
    
                # Detect cross-section start line
                if "Type RM Length L Ch R = 1" in line:
                    xs = line.split(",")[1].strip()  # Extract cross-section identifier
                    self.xs_data[xs] = {
                        "chainage": float(xs),
                        "Reach": current_reach  # Add the current reach number to the data
                    }
                    for variable in self.ice_variables:
                        self.xs_data[xs][variable] = {}
                    flag = True  # Start collecting ice parameters
    
                # Extract ice parameters
                elif flag:
                    for variable in self.ice_variables:
                        if variable in line:
                            self.xs_data[xs][variable]['val'] = line.split("=")[1].strip()
                            self.xs_data[xs][variable]['lnNum'] = count
                            if len(self.xs_data[xs]) == len(self.ice_variables):
                                flag = False  # Stop collecting once all variables are gathered

    def get_XS_manning(self):
        """
        Extracts Manning's roughness coefficient from the cross-sections.
        After detecting a '#Mann' line, it extracts the second, fifth, and eighth numbers
        from the following line as 'val_LOB', 'val_MAIN', and 'val_ROB', respectively.
    
        Args:
            None
    
        Returns:
            None
        """
        
        with open(self.geo_file, 'r') as f:
            flag = False
            for count, line in enumerate(f):
                
                if "Type RM Length L Ch R = 1" in line:
                    xs = line.split(",")[1].strip()  
                    flag = True  
    
                
                elif flag and "#Mann" in line:
                    
                    next_line = next(f).split()  
                      
                    if len(next_line) >= 8:  
                        self.xs_data[xs]['Manning'] = {
                            'val_LOB': next_line[1], 
                            'val_MAIN': next_line[4],  
                            'val_ROB': next_line[7],  
                            'lnNum': count + 2  
                        }
                    flag = False  


    def get_XS_bank_stations(self):
        """
        Extracts bank station coordinates from the cross-sections.
        After detecting a '#Mann' line, it extracts the fourth and seventh numbers
        from the following line as 'val_L' and 'val_R', respectively.
    
        Args:
            None
    
        Returns:
            None
        """
        with open(self.geo_file, 'r') as f:
            flag = False
            for count, line in enumerate(f):
                # Detect cross-section start line
                if "Type RM Length L Ch R = 1" in line:
                    xs = line.split(",")[1].strip()  
                    flag = True  #
    
                
                elif flag and "#Mann" in line:
                    
                    next_line = next(f).split()  
                    
                    if len(next_line) >= 7:  
                        self.xs_data[xs]['BankStations'] = {
                            'val_L': float(next_line[3]),  
                            'val_R': float(next_line[6]),  
                            'lnNum': count + 2  
                        }
                    flag = False  


    def get_XS_geometry(self):
        """
        Extracts geometric coordinates of cross-sections.

        Args:
            None

        Returns:
            None
        """
        with open(self.geo_file, 'r') as f:
            flag = False
            for count, line in enumerate(f):
                if "Type RM Length L Ch R = 1" in line:
                    xs = line.split(",")[1].strip()
                    self.xs_data[xs]['Geometry'] = {}
                    flag = True
                elif flag and '#Sta/Elev' in line:
                    self.xs_data[xs]['Geometry']['xy'] = []
                    flag = False


    def get_XS_main_chan_geometry(self):
        """
        Extracts the main channel geometry based on bank station data.

        Args:
            None

        Returns:
            None
        """
        for xs in self.xs_data:
            main_geometry = []
            for i in range(0, len(self.xs_data[xs]['Geometry']['xy']), 2):
                if self.xs_data[xs]['BankStations']['val_L'] <= self.xs_data[xs]['Geometry']['xy'][i] <= self.xs_data[xs]['BankStations']['val_R']:
                    main_geometry.extend(self.xs_data[xs]['Geometry']['xy'][i:i+2])
            self.xs_data[xs]['MainChannelGeometry'] = {'xy': main_geometry}

    def get_bridge_data(self):
        """
        Extracts bridge data from the geometry file.

        Args:
            None

        Returns:
            None
        """
        self.bridge = False
        self.bridge_data = {}
        with open(self.geo_file, 'r') as f:
            for count, line in enumerate(f):
                if "Type RM Length L Ch R = 3" in line:
                    self.bridge = True
                    xs = line.split(",")[1].strip()
                    self.bridge_data[str(count)] = {'chainage': float(xs)}


    def check_stochVars(self):
        """
        Checks if all required stochastic variables are provided, and raises an error if any are missing.
        """
        # List of required keys, except jam_locations which has a special case
        required_keys = ['Frontthick', 'Q', 'friction_angle', 'porosity']

        print(""" 
              \n
---------------------------------------\n
Checking format of stochastic variables\n
---------------------------------------
              """)


        # Loop through required keys and raise an error if any key is missing
        for key in required_keys:
            if key not in self.stochVars:
                raise KeyError(f"Missing required stochastic variable: '{key}'")

            # If key exists, proceed as normal
            value = self.stochVars[key]
            print(f"Variable '{key}' found with value: {value}")

        # Special handling for jam_locations, ensure it has a list provided
        if 'jam_locations' in self.stochVars:
            if not isinstance(self.stochVars['jam_locations'], list):
                raise ValueError(f"'jam_locations' must be a list of possible values.")
            print(f"Variable 'jam_locations' found with list of values: {self.stochVars['jam_locations']}")

        print('\nVariables entered correctly!')

    def generate_stochastic_values(self):
        """
        Loops through each stochastic variable in `stochVars` and generates
        a list of values according to the specified distribution type (uniform, GEV, Gumbel, or Normal).
        
        The length of the list corresponds to the number of simulations (`NSims`).
        The last value in the list represents the number of decimal places for rounding.
        """
        
        print(""" 
              \n
---------------------------------------\n
Generating stochastic variables\n
---------------------------------------
              """)
        for var, config in self.stochVars.items():
            
            # Special case for jam_locations: no distribution type, random selection from list
            if var == 'jam_locations':
                values = self.generate_random_jam_locations(config)
            else:
                # Normal case: extract distribution type, params, and decimal places
                dist_type, params, decimal_places = config
                
                if dist_type == 'uniform':
                    values = self.generate_uniform_values(params, decimal_places)
                elif dist_type == 'GEV':
                    values = self.generate_GEV_values(params[0], params[1], limit_value=params[2], decimal_places=decimal_places, plot=True)
                elif dist_type == 'Gumbel':
                    print('doing Gumbel')
                    values = self.generate_Gumbel_values(params[0], params[1], limit_value=params[2], decimal_places=decimal_places, plot=True)
                elif dist_type == 'Normal':
                    values = self.generate_normal_values(params, decimal_places)
                elif dist_type == 'Weibull':
                    values = self.generate_Weibull_values(params[0], params[1], limit_value=params[2], decimal_places=decimal_places, plot=True)
                else:
                    raise ValueError(f"Unsupported distribution type: {dist_type}")
            
            # Dynamically create the attribute for this variable
            setattr(self, var, values)
            print(f"Generated values for {var}")

    def generate_random_jam_locations(self, jam_locations_list):
        """
        Randomly selects values from the provided list of lists of jam locations for each simulation.
        
        Args:
            jam_locations_list (list): List of lists representing possible jam locations.
        
        Returns:
            list: A list of `NSims` randomly selected jam locations (each a sublist).
        """
        return [jam_locations_list[np.random.randint(len(jam_locations_list))] for _ in range(self.NSims)]


    def generate_uniform_values(self, params, decimal_places):
        """
        Generates values from a uniform distribution.
        
        Args:
            params (list): List containing the range [min, max] for the uniform distribution.
            decimal_places (int): The number of decimal places to round the generated values to.
        
        Returns:
            list: A list of `NSims` values sampled from the uniform distribution, rounded to the specified number of decimal places.
        """
        low, high = params
        values = np.random.uniform(low, high, self.NSims)
        return np.round(values, decimal_places).tolist()

    def generate_GEV_values(self, file_path, column_name, limit_value=None, decimal_places=2, plot=False, sample_size=1000):
        """
        Generates values from a GEV distribution based on historical data.
        
        Args:
            file_path (str): Path to the CSV file containing historical data.
            column_name (str): The column to extract the data from.
            limit_value (float): Optional upper limit for the values.
            decimal_places (int): The number of decimal places to round the generated values to.
            plot (bool): Whether to plot the fitted distribution against observed and sampled data.
        
        Returns:
            list: A list of `NSims` values sampled from the fitted GEV distribution, rounded to the specified number of decimal places.
        """
        self.limit_value = limit_value

        # Load and clean the data
        data = pd.read_csv(file_path)
        self.observed_values = data[column_name].dropna()

        # Fit the GEV distribution
        self.params = genextreme.fit(self.observed_values)

        # Generate random samples from the fitted GEV distribution
        self.GEV_values = genextreme.rvs(*self.params, size=self.NSims*2)

        # Apply filtering based on limit_value, if provided
        if limit_value is not None:
            self.GEV_values = self.GEV_values[self.GEV_values <= limit_value]
            
        
        # Optionally, plot the fitted distribution, observed data, and sampled data
        if plot:
            self.plot_fitted_distribution(self.observed_values, self.GEV_values, genextreme, self.params)

        self.GEV_values=np.round(self.GEV_values, decimal_places).tolist()
        self.GEV_values = random.sample(self.GEV_values, self.NSims)

        return self.GEV_values

    def generate_Gumbel_values(self, file_path, column_name, limit_value=None, decimal_places=2, plot=True, sample_size=1000):
        """
        Generates values from a Gumbel distribution based on user-supplied historical data.
        
        Args:
            file_path (str): Path to the CSV file containing historical data.
            column_name (str): The column to extract the data from.
            limit_value (float): Optional upper limit for the values.
            decimal_places (int): Number of decimal places to round the generated values.
            plot (bool): Whether to plot the fitted distribution against observed and sampled data.
        
        Returns:
            list: A list of `NSims` values sampled from the fitted Gumbel distribution, rounded to `decimal_places`.
        """
        # Load and clean the data
        data = pd.read_csv(file_path)
        historical_data = data[column_name].dropna()
    
        # Fit the Gumbel distribution (loc, scale)
        self.params = gumbel_r.fit(historical_data)
    
        # Generate random samples from the fitted Gumbel distribution
        self.Gumbel_values = gumbel_r.rvs(loc=self.params[0], scale=self.params[1], size=self.NSims*2)
    
        # Apply filtering based on limit_value, if provided
        if limit_value is not None:
            self.Gumbel_values = self.Gumbel_values[self.Gumbel_values <= limit_value]

        # Optionally, plot the fitted distribution, observed data, and sampled data
        print(plot)
        if plot:
            print('plotting')
            self.plot_fitted_distribution(historical_data, self.Gumbel_values, gumbel_r, self.params)
    
        self.Gumbel_values=np.round(self.Gumbel_values, decimal_places).tolist()
        self.Gumbel_values = random.sample(self.Gumbel_values, self.NSims)
        
        return self.Gumbel_values


    def generate_Weibull_values(self, file_path, column_name, limit_value=None, decimal_places=2, plot=False, sample_size=1000):
        """
        Generates values from a Weibull distribution based on user-supplied historical data.
        
        Args:
            file_path (str): Path to the CSV file containing historical data.
            column_name (str): The column to extract the data from.
            limit_value (float): Optional upper limit for the values.
            decimal_places (int): Number of decimal places to round the generated values.
            plot (bool): Whether to plot the fitted distribution against observed and sampled data.
        
        Returns:
            list: A list of `sample_size` values sampled from the fitted Weibull distribution, rounded to `decimal_places`.
        """
        # Load and clean the data
        data = pd.read_csv(file_path)
        historical_data = data[column_name].dropna()
        
        # Fit the Weibull distribution (shape, loc, scale)
        params = weibull_min.fit(historical_data)
        
        # Generate random samples from the fitted Weibull distribution
        weibull_values = weibull_min.rvs(*params, size=sample_size*2)
    
        # Apply filtering based on limit_value, if provided
        if limit_value is not None:
            weibull_values = weibull_values[weibull_values <= limit_value]
    
        # Optionally, plot the fitted distribution, observed data, and sampled data
        if plot:
            self.plot_fitted_distribution(historical_data, weibull_values, weibull_min, params)
    
        weibull_values = np.round(weibull_values, decimal_places).tolist()
        
        return weibull_values



    def generate_normal_values(self, params, decimal_places):
        """
        Generates values from a normal (Gaussian) distribution.
        
        Args:
            params (list): List containing [mean, standard deviation].
            decimal_places (int): The number of decimal places to round the generated values to.
        
        Returns:
            list: A list of `NSims` values sampled from the normal distribution, rounded to the specified number of decimal places.
        """
        mean, std_dev = params
        values = norm.rvs(loc=mean, scale=std_dev, size=self.NSims)
        return np.round(values, decimal_places).tolist()

    def plot_fitted_distribution(self, observed_data, sampled_data, distribution, params, bins=50, save_as_pdf=True, pdf_filename="fitted_distribution_plot.pdf"):
        """
        Plots the fitted distribution along with the observed and sampled data, and adds the distribution equation and fitted parameters.
        
        Parameters:
        - observed_data (array-like): The original observed data.
        - sampled_data (array-like): Data sampled from the fitted distribution.
        - distribution (scipy.stats distribution): The fitted distribution (e.g., scipy.stats.genextreme, scipy.stats.gumbel_r).
        - params (tuple): Parameters for the fitted distribution.
        - bins (int, optional): Number of bins for the histograms (default is 30).
        - save_as_pdf (bool, optional): If True, the plot will be saved as a PDF file. Default is False.
        - pdf_filename (str, optional): Filename for the saved PDF. Default is "fitted_distribution_plot.pdf".
        """
        
        # Create a dictionary to map the distribution to its name
        distribution_names = {
            'genextreme': 'Generalized Extreme Value',
            'gumbel_r': 'Gumbel',
            'weibull_min': 'Weibull'
        }
        
        # Get the distribution name dynamically
        dist_name = distribution.name if hasattr(distribution, 'name') else distribution_names.get(distribution.__name__, 'Unknown')
        
        # Create a figure and plot the histograms of observed and sampled data
        plt.figure(figsize=(6, 4))
        plt.hist(observed_data, bins=bins, density=True, alpha=0.6, color='black', label='Observed Data')
        plt.hist(sampled_data, bins=bins, density=True, alpha=0.6, color='red', label='Sampled Data')
        
        # Generate the x values for plotting the PDF line
        x = np.linspace(min(observed_data), max(observed_data), 1000)
        
        # Plot the fitted distribution's PDF
        pdf_fitted = distribution.pdf(x, *params)
        plt.plot(x, pdf_fitted, 'k-', lw=2, label=f'Fitted {dist_name} Distribution')
        
        # Add labels, legend, and title
        plt.xlabel('')
        plt.ylabel('Probability density')
        
        plt.legend(loc='upper right')
        
        # Display the fitted parameters on the plot
        if distribution == genextreme:
            # GEV parameters: shape, location, scale
            param_text = f"Shape (ξ): {params[0]:.3f}\nLoc (μ): {params[1]:.3f}\nScale (σ): {params[2]:.3f}"
            equation = r'$f(x; \mu, \sigma, \xi) = \frac{1}{\sigma} \left( 1 + \xi \frac{x - \mu}{\sigma} \right)^{-\frac{1}{\xi} - 1} \exp\left( - \left( 1 + \xi \frac{x - \mu}{\sigma} \right)^{-\frac{1}{\xi}} \right)$'
        elif distribution == gumbel_r:
            # Gumbel parameters: location, scale
            param_text = f"Loc (μ): {params[0]:.3f}\nScale (σ): {params[1]:.3f}"
            equation = r'$f(x; \mu, \sigma) = \frac{1}{\sigma} \exp\left( - \frac{x - \mu}{\sigma} - \exp\left(- \frac{x - \mu}{\sigma}\right) \right)$'
        elif distribution == weibull_min:
            # Weibull parameters: shape, location, scale
            param_text = f"Shape (k): {params[0]:.3f}\nLoc (μ): {params[1]:.3f}\nScale (λ): {params[2]:.3f}"
            equation = r'$f(x; k, \lambda, \mu) = \frac{k}{\lambda} \left(\frac{x - \mu}{\lambda}\right)^{k-1} \exp\left(-\left(\frac{x - \mu}{\lambda}\right)^k\right)$'
            
        # Add the parameter text and the equation to the plot
        plt.text(0.5, 0.6, param_text, transform=plt.gca().transAxes, fontsize=10, bbox=dict(facecolor='white', alpha=0.6))
        plt.text(0.5, 0.5, equation, transform=plt.gca().transAxes, fontsize=8, bbox=dict(facecolor='white', alpha=0.6))
        
        # Save the plot as a PDF if requested
        if save_as_pdf:
            plt.savefig(pdf_filename, format='pdf', bbox_inches='tight')
            print(f"Plot saved as {pdf_filename}")
        
        # Show the plot
        # plt.show()


def save_batch(stochICE_instance: StochICE_HECRAS):
    """
    Saves the current StochICE instance to a file.

    Args:
        stochICE_instance (StochICE): The current instance of StochICE to save.

    Returns:
        None
    """
    path = os.path.join(stochICE_instance.prjDir, f'Results/RIVICE/{stochICE_instance.ID}/{stochICE_instance.ID}.ice')
    with open(path, "wb") as f:
        pickle.dump(stochICE_instance, f)

def open_batch(path: str) -> StochICE_HECRAS:
    """
    Loads a saved StochICE instance from a file.

    Args:
        path (str): The file path to the saved StochICE instance.

    Returns:
        StochICE: The loaded StochICE instance.
    """
    with open(path, "rb") as f:
        return pickle.load(f)

