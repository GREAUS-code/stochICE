
import os
import shutil
import subprocess
import numpy as np
import rasterio
import glob
import re
import pandas as pd
import pickle
import statsmodels.api as sm
import matplotlib.pyplot as plt


class StochHECRASParallel:
    """
    A class to manage and execute parallel stochastic simulations
    for ice-jam flood modeling using HEC-RAS.
    """


    def __init__(self, src_folder, script_path, n, params, cleanup_after_run=False):
        """
        Constructor for the StochHECRASParallel class.
    
        Args:
            src_folder (str): Path to the source HECRAS project folder.
            script_path (str): Path to the main python script to copy and modify.
            n (int): Number of copies to create for parallel execution.
            params (dict): Configuration parameters for the simulations.
            cleanup_after_run (bool, optional): If True, cleans up generated files after execution.
        """
        import datetime
    
        self.src_folder = src_folder
        self.script_path = script_path
        self.n = n
        self.params = params
        self.cleanup_after_run = cleanup_after_run
        self.scripts_list = []
    
        # Get batch name and timestamp
        batch_name = params.get("batch_ID", "Unknown Batch")
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
        print(f"{'='*60}")
        print(f"            StochICE Parallel Simulation Launch")
        print(f"{'='*60}\n")
        print(f"Batch Name    : {batch_name}")
        print(f"Launch time   : {timestamp}")
        print(f"Processes     : {self.n}")
        print(f"Source Folder : {self.src_folder}")
        print(f"Script Path   : {self.script_path}")
        print(f"Cleanup After : {'Yes' if self.cleanup_after_run else 'No'}\n")
        print(f"Simulation Parameters:")
    
        for key, value in self.params.items():
            if isinstance(value, dict):  # Check if the value is a dictionary (e.g., stochVars)
                print(f"  {key}:")
                for sub_key, sub_value in value.items():
                    if isinstance(sub_value, list):  # If sub_value is a list, print only the first two values
                        print(f"    {sub_key}:")
                        for i, item in enumerate(sub_value[:2]):
                            print(f"      - {item}")
                        if len(sub_value) > 2:
                            pass
                    else:
                        print(f"    {sub_key}: {sub_value}")
            else:
                print(f"  {key}: {value}")
    

    
        # Initialize required directories and create scripts for simulation
        self.setup_stochice_data_directory()
        self.create_control_file()
        self.scripts_list = self.copy_folders_and_modify_script()
    
       
        print(f"\nStochICE will now run in parallel on {self.n} processes.")
        print("The console output below displays logs from the first process only:")
        
        
        # Execute simulations and process results
        if self.scripts_list:
            self.run_scripts_in_parallel()
            
            self.copy_generated_data_back()
            self.create_combined_maximum_depth_map()
            self.create_combined_frequency_map()
            self.plot_wse_envelope()

            if self.cleanup_after_run:
                self.cleanup_generated_files()
    
        # Save state and input scripts for reproducibility
        save_state(self)
        self.save_script_to_inputs()


        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
       
        print(f"\n{'='*57}")
        print(f"       Batch complete at: {timestamp}")
        print(f"{'='*57}")
        
    def create_control_file(self):
        """
        Creates a control Python script for configuring and running simulations.

        This method writes a script that initializes parameters for simulations
        based on the input `params` dictionary.

        Args:
            None

        Returns:
            None
        """
        with open(self.script_path, 'w') as file:
            file.write("import stochICE as ice\n\n")

            # Write parameters to the control file
            file.write(f"path = r'{self.params['path']}'\n")
            file.write(f"batch_ID = '{self.params['batch_ID']}'\n")
            file.write(f"ras = '{self.params['ras']}'\n")
            file.write(f"geo = '{self.params['geo']}'\n")
            file.write(f"flowFile = '{self.params['flowFile']}'\n")
            file.write(f"wse = r'{self.params['wse']}'\n\n")

            file.write(f"NSims = {self.params['NSims']}\n")
            file.write(f"ice_jam_reach = {self.params['ice_jam_reach']}\n\n")

            # Write stochastic variables
            file.write("stochVars = {\n")
            for key, value in self.params['stochVars'].items():
                file.write(f"    '{key}': {value},\n")
            file.write("}\n\n")

            # Final script setup
            file.write(
                "Chateauguay_embacle = ice.StochICE_HECRAS("
                "path, batch_ID, ras, geo, flowFile, wse, NSims, ice_jam_reach, stochVars)\n"
            )

    def setup_stochice_data_directory(self):
        """
        Sets up the directory structure for simulation data.

        If directories already exist, their contents are deleted; otherwise, new directories are created.

        Args:
            None

        Returns:
            None
        """
        self.data_path = os.path.join(self.src_folder, f"StochICE_data_{self.params['batch_ID']}")
        self.inputs_path = os.path.join(self.data_path, "Inputs")
        self.results_path = os.path.join(self.data_path, "Results")
        self.distributions_path = os.path.join(self.inputs_path, "Variable_distributions")
        self.depths_path = os.path.join(self.results_path, "Individual_depth_tifs")
        self.wse_path = os.path.join(self.results_path, "Individual_WSE_profiles")
        self.max_depth_path = os.path.join(self.results_path, "Ensemble_maximum_depth_maps")
        self.frequency_path = os.path.join(self.results_path, "Ensemble_frequency_maps")
        self.data_geo_files_path = os.path.join(self.inputs_path, "Individual_GeoFiles")
        self.data_flow_files_path = os.path.join(self.inputs_path, "Individual_FlowFiles")

        paths = [
            self.data_path,
            self.inputs_path,
            self.distributions_path,
            self.results_path,
            self.depths_path,
            self.wse_path,
            self.max_depth_path,
            self.frequency_path,
            self.data_geo_files_path,
            self.data_flow_files_path,
        ]

        for path in paths:
            if os.path.exists(path):
                # Remove existing contents
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

        print("\nStochICE data directory set up successfully.")

    def copy_folders_and_modify_script(self):
        """
        Copies the source folder for each simulation and modifies the script for parallel execution.

        This method creates `n` copies of the source folder and modifies each script to
        include the correct folder paths for the simulation.

        Args:
            None

        Returns:
            list: A list of paths to the modified scripts for each simulation.
        """
        parent_dir = os.path.dirname(self.src_folder)
        folder_name = os.path.basename(self.src_folder)
        script_dir = os.path.dirname(self.script_path)
        script_name = os.path.basename(self.script_path)
        scripts = []

        for i in range(1, self.n + 1):
            # Create a new folder for this simulation
            new_folder_name = f"{folder_name}_{i}"
            new_folder_path = os.path.join(parent_dir, new_folder_name)
            shutil.copytree(self.src_folder, new_folder_path)

            # Modify the script for the new folder
            with open(self.script_path, 'r') as script_file:
                script_content = script_file.read()
            modified_script_content = script_content.replace(
                repr(self.src_folder).strip("'"), repr(new_folder_path).strip("'")
            )

            # Save the modified script
            new_script_name = f"{os.path.splitext(script_name)[0]}_{i}.py"
            new_script_path = os.path.join(script_dir, new_script_name)
            with open(new_script_path, 'w') as new_script_file:
                new_script_file.write(modified_script_content)
            scripts.append(new_script_path)

        return scripts

    def parse_simulation_count(self, variable_name="NSims"):
        """
        Parses the simulation count from the script file by extracting the value of a specified variable.
    
        This method reads the script file specified in `self.script_path` and searches for the variable
        defined by `variable_name`. It extracts the value assigned to this variable and returns it.
    
        Args:
            variable_name (str): The name of the variable to search for in the script file. Defaults to "NSims".
    
        Returns:
            int: The value assigned to the specified variable in the script file.
    
        Raises:
            ValueError: If the specified variable is not found in the script file.
    
        Example:
            >>> count = instance.parse_simulation_count("NSims")
            >>> print(count)
            100  # Example output based on the script content.
        """
        # Read the script file to extract its content.
        with open(self.script_path, 'r') as file:
            content = file.read()
        
        # Use a regular expression to search for the variable and its value.
        match = re.search(fr'{variable_name}\s*=\s*(\d+)', content)
        if match:
            # Return the integer value of the variable.
            return int(match.group(1))
        else:
            # Raise an error if the variable is not found.
            raise ValueError(f"Could not find the '{variable_name}' variable in the script.")



    def run_scripts_in_parallel(self):
        """
        Executes simulation scripts in parallel.
    
        Real-time output of the first subprocess is streamed to the console, 
        while other subprocesses log their output to separate files.
        """
        processes = []
        log_files = []
        first_process = None
    
        for i, script in enumerate(self.scripts_list):
            log_file = open(f"process_{i + 1}.log", "w")
            log_files.append(log_file)
            
            if i == 0:
                # Capture real-time output for the first process
                first_process = subprocess.Popen(
                    ["python", "-u", script],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True
                )
            else:
                # Redirect output to a log file for other processes
                processes.append(subprocess.Popen(["python", script], stdout=log_file, stderr=subprocess.STDOUT))
    
        # Stream the first subprocess's output
        if first_process:
            try:
                for line in iter(first_process.stdout.readline, ''):
                    print(line, end='', flush=True)
            finally:
                first_process.stdout.close()
                first_process.wait()
    
        # Wait for other processes to complete
        for process in processes:
            process.wait()
    
        # Close all log files
        for log_file in log_files:
            log_file.close()


    def copy_generated_data_back(self):
        """
        Copies generated data from simulation folders back to the source folder.

        This method collects results from all parallel simulation folders
        and appends identifiers (_i) to distinguish outputs from different processes.

        Args:
            None

        Returns:
            None
        """
        data_subfolders = [
            r'Inputs\Individual_FlowFiles',
            r'Inputs\Individual_GeoFiles',
            r'Results\Individual_depth_tifs',
            r'Results\Individual_WSE_profiles',
            r'Results\Ensemble_frequency_maps',
            r'Results\Ensemble_maximum_depth_maps'
        ]

        for i in range(1, self.n + 1):
            new_folder_name = f"{os.path.basename(self.src_folder)}_{i}"
            new_folder_path = os.path.join(os.path.dirname(self.src_folder), new_folder_name)

            for subfolder in data_subfolders:
                src_subfolder_path = os.path.join(new_folder_path, f"StochICE_data_{self.params['batch_ID']}", subfolder)
                dest_subfolder_path = os.path.join(self.src_folder, f"StochICE_data_{self.params['batch_ID']}", subfolder)

                if os.path.exists(src_subfolder_path):
                    for item in os.listdir(src_subfolder_path):
                        src_item = os.path.join(src_subfolder_path, item)
                        dest_item_base = os.path.join(dest_subfolder_path, item)
                        dest_item = f"{os.path.splitext(dest_item_base)[0]}_{i}{os.path.splitext(dest_item_base)[1]}"

                        if os.path.isfile(src_item):
                            shutil.copy2(src_item, dest_item)
                        elif os.path.isdir(src_item):
                            shutil.copytree(src_item, dest_item)
                else:
                    print(f"Subfolder '{src_subfolder_path}' does not exist")

    def plot_wse_envelope(self):
        """
        Generates water surface elevation (WSE) envelope plots.

        This method groups CSV files by river number, plots the WSE profiles for
        each river, and saves the plots to the WSE_envelopes folder.

        Args:
            None

        Returns:
            None
        """
        river_files = {}

        # Collect files grouped by river number
        for root, dirs, files in os.walk(self.wse_path):
            for file in files:
                if file.endswith('.csv') and 'River' in file:
                    river_number = file.split('_')[0]  # Extract river number
                    if river_number not in river_files:
                        river_files[river_number] = []
                    river_files[river_number].append(os.path.join(root, file))

        envelope_path = os.path.join(self.results_path, "WSE_envelopes")
        os.makedirs(envelope_path, exist_ok=True)

        # Plot WSE profiles for each river
        for river_number, file_paths in river_files.items():
            plt.figure(figsize=(10, 6))
            for file_path in file_paths:
                data = pd.read_csv(file_path)
                plt.plot(data.iloc[:, 0], data.iloc[:, 1], label=os.path.basename(file_path))

            # Set plot properties
            plt.xlabel('Chainage (m)')
            plt.ylabel('Water Surface Elevation (m)')
            plt.title(f'Water Surface Elevation Profile for {river_number}')
            plt.legend()

            # Save the plot
            plot_filename = os.path.join(envelope_path, f"{river_number}_wse_envelope.png")
            plt.savefig(plot_filename)
           
            plt.close()

    def create_combined_maximum_depth_map(self):
        """
        Creates a combined maximum depth map from individual depth map files.

        This method finds and combines all maximum depth map TIFF files into
        a single raster file representing the maximum depths across all simulations.

        Args:
            None

        Returns:
            None
        """
        depth_map_files = glob.glob(os.path.join(self.max_depth_path, "maximum_depth_map_*.tif"))

        if not depth_map_files:
            print("No 'maximum_depth_map' files found.")
            return

        with rasterio.open(depth_map_files[0]) as src:
            meta = src.meta
            depth_stack = np.stack([rasterio.open(f).read(1) for f in depth_map_files])

        max_depth = np.max(depth_stack, axis=0)
        meta.update(dtype=rasterio.float32, count=1, compress='DEFLATE')

        combined_path = os.path.join(self.max_depth_path, "combined_maximum_depth_map.tif")
        with rasterio.open(combined_path, 'w', **meta) as dst:
            dst.write(max_depth.astype(rasterio.float32), 1)
        

    def create_combined_frequency_map(self):
        """
        Creates a combined flood frequency map from individual frequency map files.

        This method calculates the flood frequency across all simulations and
        generates a combined map representing the probability of flooding.

        Args:
            None

        Returns:
            None
        """
        frequency_map_files = glob.glob(os.path.join(self.frequency_path, "frequency_floodmap_*.tif"))

        if not frequency_map_files:
            print("No 'frequency_floodmap' files found.")
            return

        simulations_per_process = self.parse_simulation_count()
        total_simulations = simulations_per_process * len(frequency_map_files)

        with rasterio.open(frequency_map_files[0]) as src:
            meta = src.meta
            frequency_sum = sum(rasterio.open(f).read(1) for f in frequency_map_files)

        frequency = frequency_sum / total_simulations
        meta.update(dtype=rasterio.float32, count=1, compress='DEFLATE')

        combined_path = os.path.join(self.frequency_path, "combined_frequency_floodmap.tif")
        with rasterio.open(combined_path, 'w', **meta) as dst:
            dst.write(frequency, 1)
        

    def cleanup_generated_files(self):
        """
        Cleans up folders and scripts generated during parallel execution.

        This method deletes all copied folders and modified scripts to
        free up disk space after the simulation run.

        Args:
            None

        Returns:
            None
        """
        
        print('\n---------------------------------------------------------')
        print('                 Folder and file cleanup ')
        print('---------------------------------------------------------\n')
        
        parent_dir = os.path.dirname(self.src_folder)
        folder_name = os.path.basename(self.src_folder)
        script_dir = os.path.dirname(self.script_path)
        script_base_name = os.path.splitext(os.path.basename(self.script_path))[0]

        for i in range(1, self.n + 1):
            folder_to_remove = os.path.join(parent_dir, f"{folder_name}_{i}")
            if os.path.exists(folder_to_remove):
                shutil.rmtree(folder_to_remove)

            script_to_remove = os.path.join(script_dir, f"{script_base_name}_{i}.py")
            if os.path.exists(script_to_remove):
                os.remove(script_to_remove)
        print("Cleanup of parallel processing folders completed.")


    def probability_exceedance_curve(self, chainage, filter_string="", show_plot=False, run_regression=True):
        """
        Generates a probability of exceedance curve for a given chainage.

        This method identifies the closest chainage in the input data, calculates the exceedance
        probability for water surface elevations (WSE), and optionally performs regression analysis
        and plots the results.

        Args:
            chainage (float): The chainage value to calculate the probability exceedance curve for.
            filter_string (str): A string to filter input files by name. Default is an empty string.
            show_plot (bool): Whether to display the plot. Default is False.
            run_regression (bool): Whether to run regression analysis. Default is True.

        Returns:
            None
        """
        input_path = os.path.join(self.wse_path)

        closest_chainage = None  # Placeholder for the closest chainage

        # Find the closest chainage and set up the output directory
        for filename in os.listdir(input_path):
            if filename.endswith('.csv') and filter_string in filename:
                file_path = os.path.join(input_path, filename)
                df = pd.read_csv(file_path)
                row = df['Chainage (m)'].sub(chainage).abs().idxmin()
                closest_chainage = df.iloc[row]['Chainage (m)']
                break

        if closest_chainage is None:
            print("No files matched the filter string or chainage.")
            return

        # Prepare the output subfolder
        output_subfolder = os.path.join(
            self.data_path, "Results", "Prob_exceedance", f"{closest_chainage:.1f}"
        )
        os.makedirs(output_subfolder, exist_ok=True)

        csv_output_path = os.path.join(output_subfolder, f'{filter_string}_{closest_chainage:.1f}_PEC_data.csv')

        # Initialize lists to store data
        data_for_csv = []
        water_surface_levels = []
        probabilities = []

        # Iterate through each file to gather WSE and calculate probabilities
        for filename in os.listdir(input_path):
            if filename.endswith('.csv') and filter_string in filename:
                file_path = os.path.join(input_path, filename)
                df = pd.read_csv(file_path)
                row = df['Chainage (m)'].sub(chainage).abs().idxmin()
                wse = df.iloc[row]['wse (m)']
                water_surface_levels.append(wse)

                try:
                    parts = filename.replace('.csv', '').split('_')
                    numeric_parts = [p for p in parts if p.replace('.', '', 1).isdigit()]
                    Q, thick, phi, porosity, toe, head = map(float, numeric_parts[-6:])
                    difference = head - toe
                    data_for_csv.append([wse, Q, thick, phi, porosity, toe, head, difference])
                except ValueError:
                    print(f"Filename '{filename}' does not match the expected format.")
                    continue

        # Sort WSE for probability of exceedance calculation
        water_surface_levels.sort()
        probabilities = [1 - (count + 1) / len(water_surface_levels) for count in range(len(water_surface_levels))]

        # Save data to CSV
        if data_for_csv:
            df_output = pd.DataFrame(data_for_csv, columns=['wse', 'Q', 'thick', 'phi', 'porosity', 'toe', 'head', 'difference'])
            df_output.to_csv(csv_output_path, index=False)
            print(f"Data saved to {csv_output_path}")

            if run_regression:
                self.perform_regression_analysis(csv_output_path)

        # Plot and save the probability exceedance curve
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.scatter(probabilities, water_surface_levels, c='black', marker='o', s=0.5)
        plt.gca().invert_xaxis()
        ax.set_xscale('log')
        fig.suptitle(f'Probability of Exceedance at Chainage: {closest_chainage:.1f} m')
        ax.set_xlabel('Exceedance Probability', fontsize=10)
        ax.set_ylabel('Water Surface Elevation (m)', fontsize=10)

        output_plot_path = os.path.join(output_subfolder, f'prob_exceed_{closest_chainage:.1f}m.pdf')
        plt.savefig(output_plot_path, format='pdf', dpi=1200)
        print(f"Probability exceedance plot saved to {output_plot_path}")

        if show_plot:
            plt.show()
        else:
            plt.close()

    def perform_regression_analysis(self, file_path):
        """
        Performs regression analysis on WSE data.

        This method uses ordinary least squares (OLS) regression to analyze the relationship
        between WSE and other variables.

        Args:
            file_path (str): The path to the CSV file containing WSE data.

        Returns:
            None
        """
        df = pd.read_csv(file_path)

        y = df['wse']  # Response variable
        X = df.drop(columns=['wse'])  # Predictor variables
        X = sm.add_constant(X)  # Add intercept term

        model = sm.OLS(y, X).fit()
        regression_summary = model.summary().as_text()

        base_name = os.path.splitext(os.path.basename(file_path))[0]
        output_file_path = os.path.join(os.path.dirname(file_path), f"{base_name}_regression_results.txt")

        with open(output_file_path, 'w') as f:
            f.write(regression_summary)

        print(f"Regression results saved to {output_file_path}")

    def save_script_to_inputs(self):
        """
        Copies the script file to the Inputs directory.

        Ensures that the `script_path` is backed up in the `Inputs` folder for record-keeping.

        Args:
            None

        Returns:
            None
        """
        os.makedirs(self.inputs_path, exist_ok=True)
        script_name = os.path.basename(self.script_path)
        dest_path = os.path.join(self.inputs_path, script_name)

        try:
            shutil.copy2(self.script_path, dest_path)
            
        except Exception as e:
            print(f"Failed to copy the script to {dest_path}. Reason: {e}")


def save_state(instance):
    """
    Saves the current state of the `StochHECRASParallel` instance.

    Args:
        instance (StochHECRASParallel): The instance to save.

    Returns:
        None
    """
    directory = instance.data_path
    file_path = os.path.join(directory, f"{instance.params['batch_ID']}.ice")
    os.makedirs(directory, exist_ok=True)

    try:
        with open(file_path, 'wb') as file:
            pickle.dump(instance, file)
        
    except Exception as e:
        print(f"Failed to save state. Reason: {e}")


def load_state(file_path):
    """
    Loads the saved state of a `StochHECRASParallel` instance.

    Args:
        file_path (str): The path to the saved state file.

    Returns:
        StochHECRASParallel: The loaded instance.
    """
    with open(file_path, 'rb') as file:
        instance = pickle.load(file)
    print(f"State loaded from {file_path}")
    return instance





# import os
# import shutil
# import subprocess
# import numpy as np
# import rasterio
# import glob
# import re
# import pandas as pd
# import time
# import pickle
# import statsmodels.api as sm



# import matplotlib.pyplot as plt

# class StochHECRASParallel:

#     def __init__(self, src_folder, script_path, n, params, cleanup_after_run=False):
#         """
#         Initializes the StochHECRASParallel class with required paths and configurations.
        
#         Args:
#             src_folder (str): Path of the source folder.
#             script_path (str): Path to the main script to copy and modify.
#             n (int): Number of copies to create for parallel execution.
#             cleanup_after_run (bool): Whether to clean up copied folders and scripts after execution.
#         """
#         self.src_folder = src_folder
#         self.script_path = script_path
#         self.n = n
#         self.cleanup_after_run = cleanup_after_run
#         self.scripts_list = []
#         self.params=params
        
#         # Run initialization methods
#         # self.cleanup_previous_run()
#         self.setup_stochice_data_directory()
#         self.create_control_file()
#         self.scripts_list = self.copy_folders_and_modify_script()
        
#         if self.scripts_list:
#             self.run_scripts_in_parallel()
#             self.copy_generated_data_back()
#             self.create_combined_maximum_depth_map()
#             self.create_combined_frequency_map()
#             self.plot_wse_envelope()
            
#             if self.cleanup_after_run:
#                 self.cleanup_generated_files()

#         save_state(self)
#         self.save_script_to_inputs()
    
#     def create_control_file(self):
#         """
#         Creates a control .py file at the specified path with the provided parameters.
    
#         Args:
#             script_path (str): Path where the control file should be saved.
#             params (dict): Dictionary containing values for each variable in the control file.
#         """
#         with open(self.script_path, 'w') as file:
#             file.write("import stochICE as ice\n\n")
            
#             # Write each parameter line by line
#             file.write(f"path = r'{self.params['path']}'\n")
#             file.write(f"batch_ID = '{self.params['batch_ID']}'\n")
#             file.write(f"ras = '{self.params['ras']}'\n")
#             file.write(f"geo = '{self.params['geo']}'\n")
#             file.write(f"flowFile = '{self.params['flowFile']}'\n")
#             file.write(f"wse = r'{self.params['wse']}'\n\n")
            
#             file.write(f"NSims = {self.params['NSims']}\n")
#             file.write(f"ice_jam_reach = {self.params['ice_jam_reach']}\n\n")
            
#             # Write the stochVars dictionary with formatting
#             file.write("stochVars = {\n")
#             for key, value in self.params['stochVars'].items():
#                 if isinstance(value[1], list):
#                     # For variables that have a list (like jam_locations or distributions with parameters)
#                     file.write(f"    '{key}': {value},\n")
#                 else:
#                     file.write(f"    '{key}': {value},\n")
#             file.write("}\n\n")
            
#             # Write the final line that creates the Chateauguay_embacle object
#             file.write("Chateauguay_embacle = ice.StochICE_HECRAS(")
#             file.write("path, batch_ID, ras, geo, flowFile, wse, NSims, ice_jam_reach, stochVars)\n")
    
#     def parse_simulation_count(self, variable_name="NSims"):
#         with open(self.script_path, 'r') as file:
#             content = file.read()
#         match = re.search(fr'{variable_name}\s*=\s*(\d+)', content)
#         if match:
#             return int(match.group(1))
#         else:
#             raise ValueError(f"Could not find the '{variable_name}' variable in the script.")


#     def setup_stochice_data_directory(self):
#         """
#         Sets up directory structure to hold data related to a StochICE batch. 
#         If directories exist, their contents are erased. If not, the directories are created.
    
#         Args:
#             None
    
#         Returns:
#             None
#         """
#         self.data_path = os.path.join(self.src_folder, f"StochICE_data_{self.params['batch_ID']}")
#         self.inputs_path = os.path.join(self.data_path, "Inputs")
#         self.results_path = os.path.join(self.data_path, "Results")
        
#         self.distributions_path=os.path.join(self.inputs_path, "Variable_distributions")
        
#         self.depths_path = os.path.join(self.results_path, "Individual_depth_tifs")
#         self.wse_path = os.path.join(self.results_path, "Individual_WSE_profiles")
        
#         self.max_depth_path = os.path.join(self.results_path, "Ensemble_maximum_depth_maps")
#         self.frequency_path = os.path.join(self.results_path, "Ensemble_frequency_maps")
        
#         self.data_geo_files_path = os.path.join(self.inputs_path, "Individual_GeoFiles")
#         self.data_flow_files_path = os.path.join(self.inputs_path, "Individual_FlowFiles")
    
#         paths = [
#             self.data_path,
#             self.inputs_path,
#             self.distributions_path,
#             self.results_path,
#             self.depths_path,
#             self.wse_path,
#             self.max_depth_path,
#             self.frequency_path,
#             self.data_geo_files_path,
#             self.data_flow_files_path,
#         ]
        
#         for path in paths:
#             if os.path.exists(path):
#                 # If the directory exists, remove its contents
#                 for filename in os.listdir(path):
#                     file_path = os.path.join(path, filename)
#                     try:
#                         if os.path.isfile(file_path) or os.path.islink(file_path):
#                             os.unlink(file_path)  # Remove file or symlink
#                         elif os.path.isdir(file_path):
#                             shutil.rmtree(file_path)  # Remove directory and its contents
#                     except Exception as e:
#                         print(f"Failed to delete {file_path}. Reason: {e}")
#             else:
#                 # If the directory does not exist, create it
#                 os.makedirs(path)


#         print("StochICE data directory set up successfully.")



#     def copy_folders_and_modify_script(self):
#         parent_dir = os.path.dirname(self.src_folder)
#         folder_name = os.path.basename(self.src_folder)
#         script_dir = os.path.dirname(self.script_path)
#         script_name = os.path.basename(self.script_path)
#         scripts = []

#         for i in range(1, self.n + 1):
#             new_folder_name = f"{folder_name}_{i}"
#             new_folder_path = os.path.join(parent_dir, new_folder_name)
            
#             shutil.copytree(self.src_folder, new_folder_path)
            
#             with open(self.script_path, 'r') as script_file:
#                 script_content = script_file.read()
            
#             src_folder_raw = repr(self.src_folder).strip("'")
#             new_folder_path_raw = repr(new_folder_path).strip("'")
#             modified_script_content = script_content.replace(src_folder_raw, new_folder_path_raw)
            
#             new_script_name = f"{os.path.splitext(script_name)[0]}_{i}.py"
#             new_script_path = os.path.join(script_dir, new_script_name)
            
#             with open(new_script_path, 'w') as new_script_file:
#                 new_script_file.write(modified_script_content)
            
#             scripts.append(new_script_path)

#         return scripts


#     def run_scripts_in_parallel(self):
#         """
#         Runs the generated scripts in parallel and streams the output of the first subprocess to the console in real-time.
#         """
#         processes = []
#         first_process = None
    
#         for i, script in enumerate(self.scripts_list):
#             # For the first subprocess, capture stdout and stderr
#             if i == 0:
#                 first_process = subprocess.Popen(
#                     ["python", script],
#                     stdout=subprocess.PIPE,
#                     stderr=subprocess.STDOUT,
#                     text=True  # Ensures output is in text format, not bytes
#                 )
#             else:
#                 # For the rest, run normally
#                 processes.append(subprocess.Popen(["python", script]))
    
#         # Stream the first subprocess's output
#         if first_process:
#             for line in iter(first_process.stdout.readline, ''):
#                 print(line, end='')  # Print each line as it arrives
#             first_process.stdout.close()
#             first_process.wait()
    
#         # Wait for the other processes to finish
#         for process in processes:
#             process.wait()



#     # def run_scripts_in_parallel(self):
#     #     processes = [subprocess.Popen(["python", script]) for script in self.scripts_list]
#     #     for process in processes:
#     #         process.wait()



#     def copy_generated_data_back(self):
        
#         data_subfolders = [r'Inputs\Individual_FlowFiles', 
#                            r'Inputs\Individual_GeoFiles', 
#                            r'Results\Individual_depth_tifs', 
#                            r'Results\Individual_WSE_profiles', 
#                            r'Results\Ensemble_frequency_maps',
#                            r'Results\Ensemble_maximum_depth_maps']

#         for i in range(1, self.n + 1):
#             new_folder_name = f"{os.path.basename(self.src_folder)}_{i}"
#             new_folder_path = os.path.join(os.path.dirname(self.src_folder), new_folder_name)
            
#             for subfolder in data_subfolders:
#                 src_subfolder_path = os.path.join(new_folder_path, f"StochICE_data_{self.params['batch_ID']}", subfolder)
#                 dest_subfolder_path = os.path.join(self.src_folder, f"StochICE_data_{self.params['batch_ID']}", subfolder)
                
#                 if os.path.exists(src_subfolder_path):
#                     for item in os.listdir(src_subfolder_path):
#                         src_item = os.path.join(src_subfolder_path, item)
#                         dest_item_base = os.path.join(dest_subfolder_path, item)
                        
#                         # Append '_i' to the destination file or folder name
#                         dest_item = f"{os.path.splitext(dest_item_base)[0]}_{i}{os.path.splitext(dest_item_base)[1]}"
                                 
                        
#                         if os.path.isfile(src_item):
#                             shutil.copy2(src_item, dest_item)
#                         elif os.path.isdir(src_item):
#                             shutil.copytree(src_item, dest_item)
#                 else:
#                     print(f"Subfolder '{src_subfolder_path}' does not exist")



#     def plot_wse_envelope(self):
        
#         """
#         Finds all CSV files in a specified directory, groups them by river number, 
#         and plots the first column versus the second column for each river group.
#         Saves each plot in the WSE folder with a unique name.
#         """
        
#         # Dictionary to hold file paths grouped by river number
#         river_files = {}
    
#         # Collect files grouped by river number
#         for root, dirs, files in os.walk(self.wse_path):
#             for file in files:
#                 if file.endswith('.csv') and 'River' in file:
#                     # Extract river number from the filename (assumes format "River N_<...>.csv")
#                     river_number = file.split('_')[0]  # Get "River N" part
#                     if river_number not in river_files:
#                         river_files[river_number] = []
#                     river_files[river_number].append(os.path.join(root, file))

        
#         envelop_path = os.path.join(self.results_path, "WSE_envelopes")    
#         if not os.path.exists(envelop_path):
#             os.makedirs(envelop_path)

        
#         # Plot each river's data in a separate plot
#         for river_number, file_paths in river_files.items():
#             plt.figure(figsize=(10, 6))
            
#             # Plot each CSV file's data for the current river
#             for file_path in file_paths:
#                 data = pd.read_csv(file_path)
#                 plt.plot(data.iloc[:, 0], data.iloc[:, 1], label=os.path.basename(file_path))
            
#             # Set plot labels and title
#             plt.xlabel('Chainage (m)')
#             plt.ylabel('Water Surface Elevation (m)')
#             plt.title(f'Water Surface Elevation Profile for {river_number}')
            
#             # Save the plot with a unique filename
#             plot_filename = os.path.join(envelop_path, f"{river_number}_wse_envelope.png")
            
#             plt.savefig(plot_filename)
#             print(f"Saved plot for {river_number} as {plot_filename}")
            
#             plt.close()



#     def create_combined_maximum_depth_map(self):
        
#         # floodmaps_folder = os.path.join(self.src_folder, "FloodMaps")
#         depth_map_files = glob.glob(os.path.join(self.max_depth_path, "maximum_depth_map_*.tif"))

#         if not depth_map_files:
#             print("No 'maximum_depth_map' files found in the FloodMaps folder.")
#             return

#         with rasterio.open(depth_map_files[0]) as src:
#             meta = src.meta
#             depth_stack = np.stack([rasterio.open(f).read(1) for f in depth_map_files])

#         max_depth = np.max(depth_stack, axis=0)
#         meta.update(dtype=rasterio.float32, count=1, compress='DEFLATE')
#         combined_path = os.path.join(self.max_depth_path, "combined_maximum_depth_map.tif")

#         with rasterio.open(combined_path, 'w', **meta) as dst:
#             dst.write(max_depth.astype(rasterio.float32), 1)

#     def create_combined_frequency_map(self):
        
#         # floodmaps_folder = os.path.join(self.src_folder, "FloodMaps")
#         frequency_map_files = glob.glob(os.path.join(self.frequency_path, "frequency_floodmap_*.tif"))

#         if not frequency_map_files:
#             print("No 'frequency_floodmap' files found in the FloodMaps folder.")
#             return

#         simulations_per_process = self.parse_simulation_count()
#         total_simulations = simulations_per_process * len(frequency_map_files)

#         with rasterio.open(frequency_map_files[0]) as src:
#             meta = src.meta
#             frequency_sum = sum(rasterio.open(f).read(1) for f in frequency_map_files)

#         frequency = frequency_sum / total_simulations
#         meta.update(dtype=rasterio.float32, count=1, compress='DEFLATE')
#         combined_path = os.path.join(self.frequency_path, "combined_frequency_floodmap.tif")

#         with rasterio.open(combined_path, 'w', **meta) as dst:
#             dst.write(frequency, 1)


#     def cleanup_generated_files(self):
#         parent_dir = os.path.dirname(self.src_folder)
#         folder_name = os.path.basename(self.src_folder)
#         script_dir = os.path.dirname(self.script_path)
#         script_base_name = os.path.splitext(os.path.basename(self.script_path))[0]

#         for i in range(1, self.n + 1):
#             folder_to_remove = os.path.join(parent_dir, f"{folder_name}_{i}")
#             if os.path.exists(folder_to_remove):
#                 shutil.rmtree(folder_to_remove)

#             script_to_remove = os.path.join(script_dir, f"{script_base_name}_{i}.py")
#             if os.path.exists(script_to_remove):
#                 os.remove(script_to_remove)

#     def probability_exceedance_curve(self, chainage, filter_string="", show_plot=False, run_regression=True):
#         input_path = os.path.join(self.wse_path)
        
#         # Placeholder for the closest chainage
#         closest_chainage = None
        
#         # Iterate through files to find the closest chainage and prepare the output directory
#         for filename in os.listdir(input_path):
#             if filename.endswith('.csv') and filter_string in filename:
#                 file_path = os.path.join(input_path, filename)
                
#                 # Load the file into a DataFrame
#                 df = pd.read_csv(file_path)
                
#                 # Find the row with chainage closest to the specified chainage
#                 row = df['Chainage (m)'].sub(chainage).abs().idxmin()
#                 closest_chainage = df.iloc[row]['Chainage (m)']
#                 break  # Stop after determining the closest chainage
    
#         if closest_chainage is None:
#             print("No files matched the filter string or chainage.")
#             return
    
#         # Create a subfolder with the closest_chainage as its name
#         output_subfolder = os.path.join(self.data_path, "Results", "Prob_exceedance", f"{closest_chainage:.1f}")
#         os.makedirs(output_subfolder, exist_ok=True)
    
#         csv_output_path = os.path.join(output_subfolder, f'{filter_string}_{closest_chainage:.1f}_PEC_data.csv')
    
#         # List to collect data for the new CSV file
#         data_for_csv = []
#         water_surface_levels = []
#         probabilities = []
    
#         # Iterate through each file in the directory
#         for filename in os.listdir(input_path):
#             if filename.endswith('.csv') and filter_string in filename:
#                 file_path = os.path.join(input_path, filename)
                
#                 # Load the file into a DataFrame
#                 df = pd.read_csv(file_path)
                
#                 # Find the row with chainage closest to the specified chainage
#                 row = df['Chainage (m)'].sub(chainage).abs().idxmin()
#                 wse = df.iloc[row]['wse (m)']  # Assuming 2nd column is WSE
#                 water_surface_levels.append(wse)
                
#                 # Extract the six last numbers from the filename
#                 try:
#                     # Split the filename by underscores and filter out any non-numeric parts
#                     parts = filename.replace('.csv', '').split('_')
#                     numeric_parts = [p for p in parts if p.replace('.', '', 1).isdigit()]
#                     # Select the last six numeric parts, converting them to floats
#                     Q, thick, phi, porosity, toe, head = map(float, numeric_parts[-6:])
#                     difference = head - toe
#                     data_for_csv.append([wse, Q, thick, phi, porosity, toe, head, difference])
#                 except ValueError:
#                     print(f"Filename '{filename}' does not match the expected format.")
#                     continue
    
#         # Sort water surface levels for probability of exceedance calculation
#         water_surface_levels.sort()
        
#         # Calculate exceedance probabilities
#         for count, wse in enumerate(water_surface_levels):
#             probabilities.append(1 - (count + 1) / len(water_surface_levels))
        
#         # Save the data for CSV output
#         if data_for_csv:
#             df_output = pd.DataFrame(data_for_csv, columns=['wse', 'Q', 'thick', 'phi', 'porosity', 'toe', 'head', 'difference'])
#             df_output.to_csv(csv_output_path, index=False)
#             print(f"Data saved to {csv_output_path}")
    
#             # Run regression analysis if requested
#             if run_regression:
#                 self.perform_regression_analysis(csv_output_path)
        
#         # Generate the probability exceedance plot
#         fig, ax = plt.subplots(figsize=(6, 4))
#         ax.scatter(probabilities, water_surface_levels, c='black', marker='o', s=0.5)
        
#         # Format the plot
#         plt.gca().invert_xaxis()
#         ax.set_xscale('log')
#         fig.suptitle(f'Probability of Exceedance at Chainage: {closest_chainage:.1f} m')
#         ax.set_xlabel('Exceedance Probability', fontsize=10)
#         ax.set_ylabel('Water Surface Elevation (m)', fontsize=10)
        
#         # Save the plot
#         output_plot_path = os.path.join(output_subfolder, f'prob_exceed_{closest_chainage:.1f}m.pdf')
#         plt.savefig(output_plot_path, format='pdf', dpi=1200)
#         print(f"Probability exceedance plot saved to {output_plot_path}")
        
#         # Optionally display the plot
#         if show_plot:
#             plt.show()
#         else:
#             plt.close()


#     def perform_regression_analysis(self, file_path):
#         # Load the data
#         df = pd.read_csv(file_path)
        
#         # Define the response and predictor variables
#         y = df['wse']  # Response variable
#         X = df.drop(columns=['wse'])  # Predictor variables
#         X = sm.add_constant(X)  # Add constant term for intercept
        
#         # Perform the regression
#         model = sm.OLS(y, X).fit()
        
#         # Capture the regression summary as a string
#         regression_summary = model.summary().as_text()
        
#         # Extract the unique identifier from the file name (e.g., "River 1_1250")
#         base_name = os.path.splitext(os.path.basename(file_path))[0]
        
#         # Define the output file path with the unique identifier
#         output_file_path = os.path.join(os.path.dirname(file_path), f"{base_name}_regression_results.txt")
        
#         # Write the summary to a text file
#         with open(output_file_path, 'w') as f:
#             f.write(regression_summary)
        
#         print(f"Regression results saved to {output_file_path}")
   
        
   
#     def save_script_to_inputs(self):
#         """
#         Saves a copy of the script at `script_path` to the `Inputs` folder.
        
#         Args:
#             None
        
#         Returns:
#             None
#         """
#         # Ensure the Inputs directory exists
#         os.makedirs(self.inputs_path, exist_ok=True)
        
#         # Define the destination path for the script in the Inputs folder
#         script_name = os.path.basename(self.script_path)
#         dest_path = os.path.join(self.inputs_path, script_name)
        
#         # Copy the script to the Inputs folder
#         try:
#             shutil.copy2(self.script_path, dest_path)
#             print(f"Script copied to {dest_path}")
#         except Exception as e:
#             print(f"Failed to copy the script to {dest_path}. Reason: {e}")


# def save_state(instance):
#     """
#     Saves the current state of the StochHECRASParallel instance to a file.
    
#     Args:
#         instance (StochHECRASParallel): The instance to save.
#     """
#     # Define the directory and file path using instance attributes
#     directory = instance.data_path
#     file_path = os.path.join(directory, f"{instance.params['batch_ID']}.ice")
    
#     # Ensure the directory exists
#     os.makedirs(directory, exist_ok=True)
    
#     # Save the instance state
#     try:
#         with open(file_path, 'wb') as file:
#             pickle.dump(instance, file)
#         print(f"State saved to {file_path}")
#     except Exception as e:
#         print(f"Failed to save state. Reason: {e}")


# def load_state(file_path):
#     """
#     Loads a StochHECRASParallel instance state from a file.
    
#     Args:
#         file_path (str): Path of the saved state file.
    
#     Returns:
#         StochHECRASParallel: The loaded instance of the class.
#     """
#     with open(file_path, 'rb') as file:
#         instance = pickle.load(file)
#     print(f"State loaded from {file_path}")
#     return instance



