import os
import shutil
from datetime import datetime


import os
import sys
import shutil

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle

#import stochRIVICE
import stochHECRAS
import random
import stochHECRAS_parallel

from datetime import datetime


from scipy.stats import genextreme, gumbel_r, norm, weibull_min


class StochICE_HECRAS:
    """
    Class for managing stochastic simulations in ice-flood modeling using both HEC-RAS and RIVICE models.

    Attributes:
        prjDir (str): Directory containing the project files.
        ID (str): Batch ID for the current simulation.
        ras_file (str): Path to the HEC-RAS project file.
        geo_file (str): Path to the HEC-RAS geometry file.
        flow_file (str): Path to the HEC-RAS flow file.
        wse_map_path (str): Path to the water surface elevation map file.
        NSims (int): Number of simulations to run.
        ice_jam_reach (int): Reach ID where ice jams occur.
        stochVars (dict): Dictionary of stochastic variables and their configurations.
        xs_data (dict): Cross-section data extracted from the geometry file.
        ice_variables (list): List of ice-related variables used in the simulation.
        data_path (str): Directory for storing simulation data.
    """

    def __init__(self, prjDir: str, batch_ID: str, ras_file: str, geo_file: str,
                 flow_file: str, wse_file: str, depth_file: str, NSims: int, ice_jam_reach: int, stochvars: dict):
        """
        Initializes the StochICE_HECRAS class and sets up the simulation environment.

        Args:
            prjDir (str): Directory containing the project files.
            batch_ID (str): Batch ID for the current simulation.
            ras_file (str): Name of the HEC-RAS project file.
            geo_file (str): Name of the HEC-RAS geometry file.
            flow_file (str): Name of the HEC-RAS flow file.
            wse_file (str): Name of the water surface elevation map file.
            NSims (int): Number of simulations to run.
            ice_jam_reach (int): Reach ID where ice jams occur.
            stochvars (dict): Dictionary of stochastic variables and their configurations.
        """
        self.prjDir = prjDir
        self.ID = batch_ID
        self.ras_file = os.path.join(self.prjDir, ras_file)
        self.geo_file = os.path.join(self.prjDir, geo_file)
        self.flow_file = os.path.join(self.prjDir, flow_file)
        self.geo_file_temp = os.path.join(self.prjDir, geo_file + '.temp')
        self.wse_map_path = os.path.join(self.prjDir, "MonteCarlo", wse_file)
        self.depth_map_path = os.path.join(self.prjDir, "MonteCarlo", depth_file)

        self.NSims = NSims
        self.ice_jam_reach = ice_jam_reach
        self.stochVars = stochvars

        self.xs_data = {}
        self.ice_variables = [
            'Ice Thickness', 'Ice Mann', 'Ice Specific Gravity',
            'Ice Is Channel', 'Ice Is OB', 'Ice Friction Angle',
            'Ice Porosity', 'Ice K1', 'Ice Max Mean', 'Ice Cohesion',
            'Ice Fixed Mann'
        ]

        self.setup_stochice_data_directory()
        self.print_header()
        # self.setup_monte_carlo_dir()
        self.copy_geofile()
        self.get_XS_ice_params()
        self.get_XS_manning()
        self.get_XS_bank_stations()
        self.get_XS_geometry()
        self.get_XS_main_chan_geometry()
        self.get_bridge_data()
        self.check_stochVars()
        self.generate_stochastic_values()

        self.clr = True
        self.stochHECRAS = stochHECRAS.StochHECRAS(self)
        self.stochHECRAS.preprocess_sims()
        self.stochHECRAS.new_launch_sims()

        print("----------------------------------------------------------")
        print("            Making frequency and max depth maps")
        print("----------------------------------------------------------\n")
        self.stochHECRAS.make_frequency_flood_map()
        self.stochHECRAS.make_maximum_depth_map()
        self.stochHECRAS.make_maximum_wse_map()

    def __getstate__(self):
        """
        Customizes pickling by removing unpicklable attributes.

        Returns:
            dict: State dictionary without unpicklable attributes.
        """
        state = self.__dict__.copy()
        state['log_file'] = None
        state['original_stdout'] = None
        return state

    def __setstate__(self, state):
        """
        Restores state during unpickling.

        Args:
            state (dict): State dictionary to restore.
        """
        self.__dict__.update(state)
        self.log_file = None
        self.original_stdout = None

    def print_header(self):
        """
        Prints the header information for the current batch simulation.

        Returns:
            None
        """
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print("\n")
        print("----------------------------------------------------------")
        print("      Running stochICE flood modelling system v0.0.1      ")
        print("----------------------------------------------------------\n")
        print(f"Simulation started at: {current_time}")
        print(f"Running HEC-RAS project file: {self.ras_file}")
        print(f"HEC-RAS geometry file: {self.geo_file}")
        print(f"HEC-RAS flow file: {self.flow_file}")
        print(f"Batch ID: {self.ID}")
        print(f"Total simulations on process 1: {self.NSims}")

    def setup_stochice_data_directory(self):
        """
        Sets up directory structure to hold data related to a StochICE batch.
        If directories exist, their contents are erased; otherwise, they are created.

        Returns:
            None
        """
        self.data_path = os.path.join(self.prjDir, f"StochICE_data_{self.ID}")
        self.inputs_path = os.path.join(self.data_path, "Inputs")
        # self.outputs_path = os.path.join(self.data_path, "Outputs")
        self.results_path = os.path.join(self.data_path, "Results")

        
        self.distributions_path = os.path.join(self.inputs_path, "Variable_distributions")
        
        #depth related
        self.depth_tifs_path = os.path.join(self.results_path, "Individual_depth_tifs")
        self.max_depth_tif_path = os.path.join(self.results_path, "Ensemble_maximum_depth_maps")
        
        #water surface related
        self.wse_profiles_path = os.path.join(self.results_path, "Individual_WSE_profiles")
        self.wse_tifs_path = os.path.join(self.results_path, "Individual_WSE_tifs")
        self.max_wse_tif_path = os.path.join(self.results_path, "Ensemble_maximum_wse_maps")
        
        #frequency
        self.frequency_tif_path = os.path.join(self.results_path, "Ensemble_frequency_maps")
        
        #HECRAS files
        self.data_geo_files_path = os.path.join(self.inputs_path, "Individual_GeoFiles")
        self.data_flow_files_path = os.path.join(self.inputs_path, "Individual_FlowFiles")

        paths = [
            self.data_path, 
            self.inputs_path, 
            # self.outputs_path, 
            self.results_path, 
            self.distributions_path,
            
            self.depth_tifs_path, 
            self.max_depth_tif_path,
            self.wse_profiles_path, 
            self.wse_tifs_path, 
            self.max_wse_tif_path,
            
            self.frequency_tif_path, 
            self.data_geo_files_path, 
            self.data_flow_files_path,
        ]

        for path in paths:
            if os.path.exists(path):
                for filename in os.listdir(path):
                    file_path = os.path.join(path, filename)
                    try:
                        if os.path.isfile(file_path) or os.path.islink(file_path):
                            os.unlink(file_path)
                        elif os.path.isdir(file_path):
                            shutil.rmtree(file_path)
                    except Exception as e:
                        print(f"Failed to delete {file_path}. Reason: {e}")
            else:
                os.makedirs(path)

        
    # def setup_monte_carlo_dir(self):
    #     """
    #     Sets up directories for Monte Carlo simulations.

    #     Creates the following directories:
    #     - `SimulationTifs` for storing simulation TIF files.
    #     - `GeoFiles` for geometry files.
    #     - `FlowFiles` for flow files.

    #     Args:
    #         None

    #     Returns:
    #         None
    #     """
    #     self.MC_path = os.path.join(self.prjDir, "MonteCarlo")
    #     # self.tif_path = os.path.join(self.MC_path, "SimulationTifs")
    #     # self.geo_files_path = os.path.join(self.MC_path, "GeoFiles")
    #     # self.flow_files_path = os.path.join(self.MC_path, "FlowFiles")

    #     os.makedirs(self.MC_path, exist_ok=True)
    #     # os.makedirs(self.tif_path, exist_ok=True)
    #     # os.makedirs(self.geo_files_path, exist_ok=True)
    #     # os.makedirs(self.flow_files_path, exist_ok=True)
    #     print("Monte Carlo directory set up successfully.")

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
            print(".g0* file copied from original.\n")
        else:
            shutil.copyfile(self.geo_file, geo_file_copy)
            print("Copy of original .g0* file made.\n")

    def get_XS_ice_params(self):
        """
        Extracts ice parameters from cross-sections in the geometry file.

        Assigns each cross-section to its respective reach based on the line
        'River Reach=' (e.g., 'River Reach=River 1'). The reach number is dynamically extracted.

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
                # Check if entering a new reach
                if "River Reach=" in line:
                    # Extract the reach number
                    current_reach = int(line.split("River Reach=River")[1].split()[0].strip(","))

                # Detect cross-section start line
                if "Type RM Length L Ch R = 1" in line:
                    xs = line.split(",")[1].strip()  # Extract cross-section identifier
                    self.xs_data[xs] = {
                        "chainage": float(xs),
                        "Reach": current_reach  # Add the current reach number
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
                            # Stop collecting once all variables are gathered
                            if len(self.xs_data[xs]) == len(self.ice_variables):
                                flag = False

    def get_XS_manning(self):
        """
        Extracts Manning's roughness coefficient from cross-sections.

        After detecting a '#Mann' line, extracts the second, fifth, and eighth numbers
        from the following line as `val_LOB`, `val_MAIN`, and `val_ROB`, respectively.

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
                            'val_LOB': next_line[1],  # Left Overbank
                            'val_MAIN': next_line[4],  # Main Channel
                            'val_ROB': next_line[7],  # Right Overbank
                            'lnNum': count + 2
                        }
                    flag = False

    def get_XS_bank_stations(self):
        """
        Extracts bank station coordinates from cross-sections.

        After detecting a '#Mann' line, extracts the fourth and seventh numbers
        from the following line as `val_L` and `val_R`, respectively.

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

                    if len(next_line) >= 7:
                        self.xs_data[xs]['BankStations'] = {
                            'val_L': float(next_line[3]),  # Left Bank
                            'val_R': float(next_line[6]),  # Right Bank
                            'lnNum': count + 2
                        }
                    flag = False

    def get_XS_geometry(self):
        """
        Extracts geometric coordinates of cross-sections.

        The geometry is stored in the 'xy' key for each cross-section.

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

        Iterates over the geometric coordinates of each cross-section and selects
        points within the bank station boundaries.

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

        Bridge cross-sections are identified by the line "Type RM Length L Ch R = 3".
        The extracted data is stored in `self.bridge_data` with the line number as the key
        and chainage as the value.

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
        Checks if all required stochastic variables are provided.

        This method ensures that the required stochastic variables are present in `self.stochVars`.
        If any required variable is missing, a KeyError is raised. Additionally, it validates
        that 'jam_locations', if provided, is a list.

        Args:
            None

        Returns:
            None

        Raises:
            KeyError: If a required stochastic variable is missing.
            ValueError: If 'jam_locations' is not a list.
        """
        # List of required keys, except 'jam_locations', which has a special case
        required_keys = ['Frontthick', 'Q', 'friction_angle', 'porosity','ice_cover_n']

        print('-----------------------------------------------------------')
        print('          Checking format of stochastic variables')
        print('-----------------------------------------------------------\n')

        # Validate required keys
        for key in required_keys:
            if key not in self.stochVars:
                raise KeyError(f"Missing required stochastic variable: '{key}'")

            # Log found variable
            value = self.stochVars[key]
            print(f"Variable '{key}' found with value: {value}")

        # Special handling for 'jam_locations'
        if 'jam_locations' in self.stochVars:
            if not isinstance(self.stochVars['jam_locations'], list):
                raise ValueError("'jam_locations' must be a list of possible values.")
            print(f"Variable 'jam_locations' found with list of values: {self.stochVars['jam_locations']}")

        print('\nVariables entered correctly!\n')

    def generate_stochastic_values(self):
        """
        Generates stochastic values for each variable in `self.stochVars`.

        This method iterates through the stochastic variables and generates a list of values
        based on the specified distribution type (e.g., uniform, GEV, Gumbel, Normal, Weibull).
        The generated values are rounded to the specified number of decimal places and stored
        as attributes of the object.

        Args:
            None

        Returns:
            None

        Raises:
            ValueError: If an unsupported distribution type is encountered.
        """
        print('-----------------------------------------------------------')
        print('               Generating stochastic variables')
        print('-----------------------------------------------------------\n')

        for var, config in self.stochVars.items():
            # Handle 'jam_locations' separately
            if var == 'jam_locations':
                values = self.generate_random_jam_locations(config)
            else:
                # Extract distribution type, parameters, and decimal places
                dist_type, params, decimal_places = config

                if dist_type == 'uniform':
                    values = self.generate_uniform_values(var, params, decimal_places)
                elif dist_type == 'GEV':
                    values = self.generate_GEV_values(var, params[0], params[1], limit_value=params[2],
                                                      decimal_places=decimal_places, plot=True)
                elif dist_type == 'Gumbel':
                    print('Generating Gumbel values')
                    values = self.generate_Gumbel_values(var, params[0], params[1], limit_value=params[2],
                                                         decimal_places=decimal_places, plot=True)
                elif dist_type == 'Normal':
                    values = self.generate_normal_values(params, decimal_places)
                elif dist_type == 'Weibull':
                    values = self.generate_Weibull_values(var, params[0], params[1], limit_value=params[2],
                                                          decimal_places=decimal_places, plot=True)
                else:
                    raise ValueError(f"Unsupported distribution type: {dist_type}")

            # Dynamically create the attribute for the variable
            setattr(self, var, values)
            print(f"Generated values for {var}")

    def generate_random_jam_locations(self, jam_locations_list):
        """
        Randomly selects values from a list of jam location possibilities for each simulation.

        Args:
            jam_locations_list (list): List of lists representing possible jam locations.

        Returns:
            list: A list of `NSims` randomly selected jam locations (each as a sublist).
        """
        return [jam_locations_list[np.random.randint(len(jam_locations_list))] for _ in range(self.NSims)]

    def generate_uniform_values(self, var, params, decimal_places):
        """
        Generates values from a uniform distribution.

        Args:
            var (str): The name of the stochastic variable.
            params (list): List containing the range [min, max] for the uniform distribution.
            decimal_places (int): The number of decimal places to round the generated values to.

        Returns:
            list: A list of `NSims` values sampled from the uniform distribution,
                  rounded to the specified number of decimal places.
        """
        low, high = params
        values = np.random.uniform(low, high, self.NSims)
        return np.round(values, decimal_places).tolist()


    def generate_GEV_values(self, var, file_path, column_name, limit_value=None, decimal_places=2, plot=False, sample_size=1000):
        """
        Generates values from a Generalized Extreme Value (GEV) distribution based on historical data.

        Args:
            var (str): The name of the stochastic variable.
            file_path (str): Path to the CSV file containing historical data.
            column_name (str): The column to extract the data from.
            limit_value (float, optional): Optional upper limit for the generated values. Defaults to None.
            decimal_places (int, optional): The number of decimal places to round the generated values to. Defaults to 2.
            plot (bool, optional): Whether to plot the fitted distribution against observed and sampled data. Defaults to False.
            sample_size (int, optional): Number of samples to generate for plotting purposes. Defaults to 1000.

        Returns:
            list: A list of `NSims` values sampled from the fitted GEV distribution, rounded to the specified decimal places.
        """
        self.limit_value = limit_value

        # Load and clean the data
        data = pd.read_csv(file_path)
        self.observed_values = data[column_name].dropna()

        # Fit the GEV distribution
        self.params = genextreme.fit(self.observed_values)

        # Generate random samples
        self.GEV_values = genextreme.rvs(*self.params, size=self.NSims * 2)

        # Apply filtering if limit_value is specified
        if limit_value is not None:
            self.GEV_values = self.GEV_values[self.GEV_values <= limit_value]

        # Optionally plot the fitted distribution
        if plot:
            self.plot_fitted_distribution(var, self.observed_values, self.GEV_values, genextreme, self.params)

        # Round and downsample the values
        self.GEV_values = np.round(self.GEV_values, decimal_places).tolist()
        self.GEV_values = random.sample(self.GEV_values, self.NSims)

        return self.GEV_values

    def generate_Gumbel_values(self, var, file_path, column_name, limit_value=None, decimal_places=2, plot=True, sample_size=1000):
        """
        Generates values from a Gumbel distribution based on historical data.

        Args:
            var (str): The name of the stochastic variable.
            file_path (str): Path to the CSV file containing historical data.
            column_name (str): The column to extract the data from.
            limit_value (float, optional): Optional upper limit for the generated values. Defaults to None.
            decimal_places (int, optional): The number of decimal places to round the generated values to. Defaults to 2.
            plot (bool, optional): Whether to plot the fitted distribution against observed and sampled data. Defaults to True.
            sample_size (int, optional): Number of samples to generate for plotting purposes. Defaults to 1000.

        Returns:
            list: A list of `NSims` values sampled from the fitted Gumbel distribution, rounded to the specified decimal places.
        """
        # Load and clean the data
        data = pd.read_csv(file_path)
        historical_data = data[column_name].dropna()

        # Fit the Gumbel distribution
        self.params = gumbel_r.fit(historical_data)

        # Generate random samples
        self.Gumbel_values = gumbel_r.rvs(loc=self.params[0], scale=self.params[1], size=self.NSims * 2)

        # Apply filtering if limit_value is specified
        if limit_value is not None:
            self.Gumbel_values = self.Gumbel_values[self.Gumbel_values <= limit_value]

        # Optionally plot the fitted distribution
        if plot:
            self.plot_fitted_distribution(var, historical_data, self.Gumbel_values, gumbel_r, self.params)

        # Round and downsample the values
        self.Gumbel_values = np.round(self.Gumbel_values, decimal_places).tolist()
        self.Gumbel_values = random.sample(self.Gumbel_values, self.NSims)

        return self.Gumbel_values

    def generate_Weibull_values(self, var, file_path, column_name, limit_value=None, decimal_places=2, plot=False, sample_size=1000):
        """
        Generates values from a Weibull distribution based on historical data.

        Args:
            var (str): The name of the stochastic variable.
            file_path (str): Path to the CSV file containing historical data.
            column_name (str): The column to extract the data from.
            limit_value (float, optional): Optional upper limit for the generated values. Defaults to None.
            decimal_places (int, optional): The number of decimal places to round the generated values to. Defaults to 2.
            plot (bool, optional): Whether to plot the fitted distribution against observed and sampled data. Defaults to False.
            sample_size (int, optional): Number of samples to generate for plotting purposes. Defaults to 1000.

        Returns:
            list: A list of `sample_size` values sampled from the fitted Weibull distribution, rounded to the specified decimal places.
        """
        # Load and clean the data
        data = pd.read_csv(file_path)
        historical_data = data[column_name].dropna()

        # Fit the Weibull distribution
        params = weibull_min.fit(historical_data)

        # Generate random samples
        weibull_values = weibull_min.rvs(*params, size=sample_size * 2)

        # Apply filtering if limit_value is specified
        if limit_value is not None:
            weibull_values = weibull_values[weibull_values <= limit_value]

        # Optionally plot the fitted distribution
        if plot:
            self.plot_fitted_distribution(var, historical_data, weibull_values, weibull_min, params)

        # Round the values
        weibull_values = np.round(weibull_values, decimal_places).tolist()

        return weibull_values

    def generate_normal_values(self, params, decimal_places):
        """
        Generates values from a normal (Gaussian) distribution.

        Args:
            params (list): List containing [mean, standard deviation].
            decimal_places (int): The number of decimal places to round the generated values to.

        Returns:
            list: A list of `NSims` values sampled from the normal distribution, rounded to the specified decimal places.
        """
        mean, std_dev = params
        values = norm.rvs(loc=mean, scale=std_dev, size=self.NSims)
        return np.round(values, decimal_places).tolist()


    def plot_fitted_distribution(self, var, observed_data, sampled_data, distribution, params, bins=50, save_as_pdf=True, pdf_filename="fitted_distribution_plot.pdf"):
        """
        Plots the fitted distribution along with the observed and sampled data, and displays the distribution equation with fitted parameters.
    
        Args:
            var (str): The variable name being plotted.
            observed_data (array-like): The original observed data.
            sampled_data (array-like): Data sampled from the fitted distribution.
            distribution (scipy.stats distribution): The fitted distribution (e.g., scipy.stats.genextreme, scipy.stats.gumbel_r).
            params (tuple): Parameters for the fitted distribution.
            bins (int, optional): Number of bins for the histograms (default is 50).
            save_as_pdf (bool, optional): If True, saves the plot as a PDF. Default is True.
            pdf_filename (str, optional): The filename for the saved PDF. Default is "fitted_distribution_plot.pdf".
    
        Returns:
            None
        """
        # Map distributions to human-readable names
        distribution_names = {
            'genextreme': 'Generalized Extreme Value',
            'gumbel_r': 'Gumbel',
            'weibull_min': 'Weibull'
        }
    
        # Get the distribution name dynamically
        dist_name = distribution.name if hasattr(distribution, 'name') else distribution_names.get(distribution.__name__, 'Unknown')
    
        # Create a histogram for observed and sampled data
        plt.figure(figsize=(6, 4))
        plt.hist(observed_data, bins=bins, density=True, alpha=0.6, color='black', label='Observed Data')
        plt.hist(sampled_data, bins=bins, density=True, alpha=0.6, color='red', label='Sampled Data')
    
        # Generate x values for the fitted PDF
        x = np.linspace(min(observed_data), max(observed_data), 1000)
        pdf_fitted = distribution.pdf(x, *params)
        plt.plot(x, pdf_fitted, 'k-', lw=2, label=f'Fitted {dist_name} Distribution')
    
        # Add labels and legend
        plt.xlabel('Value')
        plt.ylabel('Probability Density')
        plt.legend(loc='upper right')
    
        # Display the fitted parameters and equation on the plot
        if distribution == genextreme:
            param_text = f"Shape (ξ): {params[0]:.3f}\nLoc (μ): {params[1]:.3f}\nScale (σ): {params[2]:.3f}"
            equation = r'$f(x; \mu, \sigma, \xi) = \frac{1}{\sigma} \left( 1 + \xi \frac{x - \mu}{\sigma} \right)^{-\frac{1}{\xi} - 1} \exp\left( - \left( 1 + \xi \frac{x - \mu}{\sigma} \right)^{-\frac{1}{\xi}} \right)$'
        elif distribution == gumbel_r:
            param_text = f"Loc (μ): {params[0]:.3f}\nScale (σ): {params[1]:.3f}"
            equation = r'$f(x; \mu, \sigma) = \frac{1}{\sigma} \exp\left( - \frac{x - \mu}{\sigma} - \exp\left(- \frac{x - \mu}{\sigma}\right) \right)$'
        elif distribution == weibull_min:
            param_text = f"Shape (k): {params[0]:.3f}\nLoc (μ): {params[1]:.3f}\nScale (λ): {params[2]:.3f}"
            equation = r'$f(x; k, \lambda, \mu) = \frac{k}{\lambda} \left(\frac{x - \mu}{\lambda}\right)^{k-1} \exp\left(-\left(\frac{x - \mu}{\lambda}\right)^k\right)$'
    
        plt.text(0.5, 0.6, param_text, transform=plt.gca().transAxes, fontsize=10, bbox=dict(facecolor='white', alpha=0.6))
        plt.text(0.5, 0.5, equation, transform=plt.gca().transAxes, fontsize=8, bbox=dict(facecolor='white', alpha=0.6))
    
        # Save plot as a PDF if requested
        if save_as_pdf:
            save_path = os.path.join(self.distributions_path, f"{var}_distribution.pdf")
            plt.savefig(save_path, format='pdf', bbox_inches='tight')
            # print(f"Plot saved as {save_path}")
    
        # Uncomment to display the plot during runtime
        # plt.show()
        plt.close()
    
    def start_logging(self, log_filename="simulation_log.txt"):
        """
        Starts logging console output to a specified file.
    
        Args:
            log_filename (str): The filename for logging console output.
    
        Returns:
            None
        """
        self.log_file = open(os.path.join(self.outputs_path, log_filename), "w")
        self.original_stdout = sys.stdout  # Save the original stdout
        sys.stdout = self.log_file  # Redirect stdout to the log file
        print(f"Logging started at {datetime.now()}")  # Start log with a timestamp
    
    def stop_logging(self):
        """
        Stops logging and restores console output to the original state.
    
        Args:
            None
    
        Returns:
            None
        """
        print(f"Logging stopped at {datetime.now()}")  # End log with a timestamp
        sys.stdout = self.original_stdout  # Restore original stdout
        self.log_file.close()  # Close the log file


def save_state(instance):

    path = os.path.join(instance.outputs_path, f'{instance.ID}.ice')
    with open(path, "wb") as f:
        pickle.dump(instance, f)

def load_state(file_path):

    with open(file_path, 'rb') as file:
        instance = pickle.load(file)
    print(f"State loaded from {file_path}")
    return instance
