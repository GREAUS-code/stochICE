
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
import datetime
import sys
import select
import time

import threading


from rasterio.features import shapes
import geopandas as gpd
from shapely.geometry import shape



from shapely.ops import unary_union


class StochHECRASParallel:
    """
    A class to manage and execute parallel stochastic simulations
    for ice-jam flood modeling using HEC-RAS.
    """

    def __init__(self, src_folder, script_path, n, params, cleanup_after_run=False, seed=[]):
        """
        Constructor for the StochHECRASParallel class.
    
        Args:
            src_folder (str): Path to the source HECRAS project folder.
            script_path (str): Path to the main python script to copy and modify.
            n (int): Number of copies to create for parallel execution.
            params (dict): Configuration parameters for the simulations.
            cleanup_after_run (bool, optional): If True, cleans up generated files after execution.
        """
        
        self.src_folder = src_folder
        self.script_path = script_path
        self.n = n
        self.params = params
        self.cleanup_after_run = cleanup_after_run
        self.scripts_list = []
    
        #location of list of parameters to resuse
        if seed:
            self.seed=seed
            print(self.seed[0])
            print(self.seed[1])
    
        # Get batch name and timestamp
        batch_name = params.get("batch_ID", "Unknown Batch")
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
        print(f"\n{'='*60}")
        print(f"            StochICE Parallel Simulation Launch")
        print(f"{'='*60}\n")
        print(f"Batch Name    : {batch_name}")
        print(f"Launch time   : {timestamp}")
        print(f"Processes     : {self.n}")
        print(f"Source Folder : {self.src_folder}")
        print(f"Script Path   : {self.script_path}")
        print(f"Cleanup After : {'Yes' if self.cleanup_after_run else 'No'}\n")
        print(f"Simulation Parameters:\n")
    
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
        
        self.create_control_file()
        self.setup_stochice_data_directory()
        self.scripts_list = self.copy_folders_and_modify_script()
        
        print(f"\nStochICE will now run in parallel on {self.n} processes.")
        print("The console output below displays logs from the first process only:")
        
        # Execute simulations and process results
        if self.scripts_list:
            
            # self.setup_stochice_data_directory()
            self.run_scripts_in_parallel()
            self.copy_generated_data_back()
            self.create_combined_maximum_depth_map()
            self.create_combined_frequency_map()
            self.create_combined_maximum_wse_map()
            # self.make_matrice_intensite()
            self.process_objectives_de_protection()
            self.plot_wse_envelope()

            if self.cleanup_after_run:
                self.cleanup_generated_files()
    
        # Save state and input scripts for reproducibility
        # save_state(self)
        # self.save_script_to_inputs()

        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
       
        print(f"\n{'='*57}")
        print(f"       Batch complete at: {timestamp}")
        print(f"{'='*57}\n")
        
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
            file.write(f"wse = '{self.params['wse']}'\n\n")
            file.write(f"depth = '{self.params['depth']}'\n\n")

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
                "path, batch_ID, ras, geo, flowFile, wse, depth, NSims, ice_jam_reach, stochVars)\n"
            )

    def setup_stochice_data_directory(self):
        """
        Sets up the directory structure for simulation data.
    
        If the directory already exists, prompts the user whether to overwrite it or exit.
    
        Args:
            None
    
        Returns:
            None
        """
        self.data_path = os.path.join(self.src_folder, f"StochICE_data_{self.params['batch_ID']}")
        self.inputs_path = os.path.join(self.data_path, "Inputs")
        self.results_path = os.path.join(self.data_path, "Results")
        self.logs_path = os.path.join(self.data_path, "Logs")
        self.distributions_path = os.path.join(self.inputs_path, "Variable_distributions")
        self.depth_tifs_path = os.path.join(self.results_path, "Individual_depth_tifs")
        self.wse_tifs_path = os.path.join(self.results_path, "Individual_wse_tifs")
        self.wse_profiles_path = os.path.join(self.results_path, "Individual_WSE_profiles")
        self.max_wse_path = os.path.join(self.results_path, "Ensemble_maximum_wse_maps")
        self.max_depth_path = os.path.join(self.results_path, "Ensemble_maximum_depth_maps")
        self.frequency_tif_path = os.path.join(self.results_path, "Ensemble_frequency_maps")
        self.obj_prot_path = os.path.join(self.results_path, "Objectives_de_protection")
        self.data_geo_files_path = os.path.join(self.inputs_path, "Individual_GeoFiles")
        self.data_flow_files_path = os.path.join(self.inputs_path, "Individual_FlowFiles")
    
        paths = [
            self.data_path,
            self.inputs_path,
            self.distributions_path,
            self.results_path,
            self.logs_path,
            self.depth_tifs_path,
            self.wse_tifs_path,
            self.wse_profiles_path,
            self.max_depth_path,
            self.max_wse_path,
            self.frequency_tif_path,
            self.obj_prot_path,
            self.data_geo_files_path,
            self.data_flow_files_path,
        ]
    
        if os.path.exists(self.data_path):
            while True:

                print(f"\n{'='*57}")
                print(f"                   User input required")
                print(f"{'='*57}")
                
                user_input = input(f"\nThe directory {self.data_path} already exists. Overwrite? (y/N): ").strip().lower()
                if user_input == 'y':
                    shutil.rmtree(self.data_path)
                    break
                elif user_input == 'n' or user_input == '':
                    print("Operation cancelled. Exiting StochICE.")
                    exit()
                else:
                    print("Invalid input. Please enter 'y' to overwrite or 'n' to exit.")
    
        for path in paths:
            os.makedirs(path, exist_ok=True)
        
        print("\nStochICE data directory set up successfully.")


    def copy_folders_and_modify_script(self):
        """
        Copies the source folder for each simulation and copies the launch script for parallel execution.
        Also splits a CSV file containing simulation parameters into `n` smaller files,
        ensuring each split file contains the header.
        
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
    
        src_folder_abs = os.path.abspath(self.src_folder)
        script_path_abs = os.path.abspath(self.script_path)
    
        existing_found = False
        for i in range(1, self.n + 1):
            new_folder_name = f"{folder_name}_{i}"
            new_folder_path = os.path.join(parent_dir, new_folder_name)
            new_script_name = f"{os.path.splitext(script_name)[0]}_{i}.py"
            new_script_path = os.path.join(script_dir, new_script_name)
    
            if os.path.exists(new_folder_path) or os.path.exists(new_script_path):
                existing_found = True
                break
    
        if existing_found:
            print(f"\n{'='*57}")
            print(f"                   User input required")
            print(f"{'='*57}")
            user_input = input("\nSome parallel folders or scripts already exist. Permission to delete them before proceeding? (y/N): ").strip().lower()
            if user_input != 'y':
                print("\nOperation cancelled. Exiting StochICE.\n")
                exit()
    
        for i in range(1, self.n + 1):
            new_folder_name = f"{folder_name}_{i}"
            new_folder_path = os.path.join(parent_dir, new_folder_name)
            new_script_name = f"{os.path.splitext(script_name)[0]}_{i}.py"
            new_script_path = os.path.join(script_dir, new_script_name)
    
            if os.path.exists(new_folder_path):
                shutil.rmtree(new_folder_path)
            if os.path.exists(new_script_path):
                os.remove(new_script_path)
    
        print(f"\n{'='*57}")
        print(f"             Creating parallel folders")
        print(f"{'='*57}\n")
    
        for i in range(1, self.n + 1):
            new_folder_name = f"{folder_name}_{i}"
            new_folder_path = os.path.join(parent_dir, new_folder_name)
    
            def ignore_stochice_data(dir, files):
                return [f for f in files if os.path.isdir(os.path.join(dir, f)) and f.startswith('StochICE_data')]
    
            shutil.copytree(src_folder_abs, new_folder_path, ignore=ignore_stochice_data)
    
            with open(script_path_abs, 'r') as script_file:
                script_content = script_file.read()
    
            modified_script_content = script_content.replace(src_folder_abs, new_folder_path)
    
            new_script_name = f"{os.path.splitext(script_name)[0]}_{i}.py"
            new_script_path = os.path.join(script_dir, new_script_name)
            with open(new_script_path, 'w') as new_script_file:
                new_script_file.write(modified_script_content)
            scripts.append(new_script_path)


        if hasattr(self, 'seed') and isinstance(self.seed, list) and len(self.seed) > 0:
            csv_path = self.seed[0]
            csv_filename = os.path.basename(csv_path)
            if os.path.exists(csv_path):
                df = pd.read_csv(csv_path)
                num_rows = len(df)
                
                if num_rows % self.n != 0:
                    valid_splits = [i for i in range(1, num_rows + 1) if num_rows % i == 0]
                    print(f"Error: The number of data lines ({num_rows}) in '{csv_filename}' cannot be evenly divided into {self.n} processes.")
                    print(f"Consider using one of the following values for 'n' to allow an even split: {valid_splits}.")
                    exit()
                
                if isinstance(self.params, dict) and 'NSims' in self.params:
                    expected_rows = self.params['NSims'] * self.n
                    if expected_rows != num_rows:
                        correct_NSims = num_rows // self.n
                        print(f"Error: The total number of data lines (i.e. {num_rows}) in '{csv_filename}' does not equal NSims*n_procs (i.e. {expected_rows}).")
                        print(f"Consider setting 'NSims' to {correct_NSims} to ensure consistency.")
                        exit()
        
                chunk_size = num_rows // self.n
                start = 0
        
                for i in range(1, self.n + 1):
                    end = start + chunk_size
                    df_chunk = df.iloc[start:end]
                    start = end
        
                    new_folder_name = f"{folder_name}_{i}"
                    new_folder_path = os.path.join(parent_dir, new_folder_name)
                    chunk_csv_path = os.path.join(new_folder_path, "seed_parameters.csv")
                    
                    with open(chunk_csv_path, 'w', newline='') as f:
                        for item in self.seed[1]:
                            f.write(f"{item}\n")
                        df_chunk.to_csv(f, index=False, mode='a', header=True)

            else:
                print(f"Error: Seed CSV file '{csv_filename}' not found. Exiting StochICE.")
                exit()
        
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
    
        - The first subprocess streams output to the console in real-time.
        - Other subprocesses log their output to separate files.
        """
        processes = []
        log_files = []
        first_process = None
    
        # Ensure logs directory exists
        os.makedirs(self.logs_path, exist_ok=True)
    
        # Ensure stdout is unbuffered for real-time console output
        sys.stdout.reconfigure(line_buffering=True)
    
        for i, script in enumerate(self.scripts_list):
            log_file_path = os.path.join(self.logs_path, f"process_{i + 1}.log")
            log_file = open(log_file_path, "w")
            log_files.append(log_file)
    
            if i == 0:
                # Start the first process with real-time output
                first_process = subprocess.Popen(
                    ["python", "-u", script],  # -u ensures unbuffered output
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    bufsize=1
                )
    
            else:
                # Other processes log output to files
                processes.append(subprocess.Popen(["python", "-u", script], stdout=log_file, stderr=subprocess.STDOUT))
    
        # Process first script's output in real-time
        if first_process:
            while True:
                output = first_process.stdout.readline()
                if output:
                    print(output, end='', flush=True)  # Print immediately
                if first_process.poll() is not None and not output:
                    break  # Exit loop if process has finished
    
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
        Copies generated data from simulation folders back to the source folder and
        vertically stacks the contents of `sim_parameters.csv` files from each simulation folder.
    
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
            r'Results\Individual_WSE_tifs',
            r'Results\Ensemble_frequency_maps',
            r'Results\Ensemble_maximum_depth_maps',
            r'Results\Ensemble_maximum_wse_maps'
        ]
        
        all_dataframes = []
    
        for i in range(1, self.n + 1):
            new_folder_name = f"{os.path.basename(self.src_folder)}_{i}"
            new_folder_path = os.path.join(os.path.dirname(self.src_folder), new_folder_name)
            
            # Process `sim_parameters.csv` if it exists
            csv_path = os.path.join(new_folder_path, "sim_parameters.csv")
            if os.path.exists(csv_path):
                df = pd.read_csv(csv_path)
                df['simulation_id'] = i  # Add identifier column
                all_dataframes.append(df)
            
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
        
        # Stack and save the consolidated `sim_parameters.csv`
        if all_dataframes:
            combined_df = pd.concat(all_dataframes, ignore_index=True)
            combined_csv_path = os.path.join(self.logs_path, "Batch_sim_parameters.csv")
            combined_df.to_csv(combined_csv_path, index=False)
            print(f"Batch simulation parameters saved to {combined_csv_path}")


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
        for root, dirs, files in os.walk(self.wse_profiles_path):
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
            print("Warning: No 'maximum_depth_map' files found. Path too long?")
            return

        with rasterio.open(depth_map_files[0]) as src:
            meta = src.meta
            depth_stack = np.stack([rasterio.open(f).read(1) for f in depth_map_files])

        max_depth = np.max(depth_stack, axis=0)
        meta.update(dtype=rasterio.float32, count=1, compress='DEFLATE')

        combined_path = os.path.join(self.max_depth_path, "combined_maximum_depth_map.tif")
        with rasterio.open(combined_path, 'w', **meta) as dst:
            dst.write(max_depth.astype(rasterio.float32), 1)
        

    def create_combined_maximum_wse_map(self):
        """
        Creates a combined maximum wse map from individual wse map files.

        This method finds and combines all maximum wse map TIFF files into
        a single raster file representing the maximum wse across all simulations.

        Args:
            None

        Returns:
            None
        """
        wse_map_files = glob.glob(os.path.join(self.max_wse_path, "maximum_wse_map_*.tif"))

        if not wse_map_files:
            print("No 'maximum_wse_map' files found.")
            return

        with rasterio.open(wse_map_files[0]) as src:
            meta = src.meta
            wse_stack = np.stack([rasterio.open(f).read(1) for f in wse_map_files])

        max_wse = np.max(wse_stack, axis=0)
        meta.update(dtype=rasterio.float32, count=1, compress='DEFLATE')

        combined_path = os.path.join(self.max_wse_path, "combined_maximum_wse_map.tif")
        with rasterio.open(combined_path, 'w', **meta) as dst:
            dst.write(max_wse.astype(rasterio.float32), 1)
        

    def create_combined_frequency_map(self):
        """
        Creates a combined flood frequency map from individual frequency map files.
    
        This method calculates the global flood frequency across all simulations and
        generates a combined map with values representing the percentage of flooding (0-100).
    
        Args:
            None
    
        Returns:
            None
        """
        import numpy as np
    
        # Get list of frequency map files
        frequency_map_files = glob.glob(os.path.join(self.frequency_tif_path, "frequency_floodmap_*.tif"))
        if not frequency_map_files:
            print("No 'frequency_floodmap' files found.")
            return
    
        # Number of simulations per map (from your parse method)
        simulations_per_process = self.parse_simulation_count()
        num_maps = len(frequency_map_files)
        total_simulations = simulations_per_process * num_maps
    
        # Read and process the maps
        with rasterio.open(frequency_map_files[0]) as src:
            meta = src.meta
            # Initialize array to accumulate the number of flooded instances
            flooded_count = np.zeros(src.shape, dtype=np.float32)
    
        # For each map, convert percentage to flooded count and accumulate
        for f in frequency_map_files:
            with rasterio.open(f) as src:
                frequency_data = src.read(1)  # Values are 0-100 (percentage)
                # Convert percentage to number of flooded simulations for this map
                flooded_in_map = (frequency_data / 100.0) * simulations_per_process
                flooded_count += flooded_in_map
    
        # Calculate global frequency as a percentage
        combined_frequency = (flooded_count / total_simulations) * 100.0
        combined_frequency = np.clip(combined_frequency, 0, 100)  # Ensure range 0-100
    
        # Update metadata and write output
        meta.update(dtype=rasterio.float32, count=1, compress='DEFLATE')
        combined_path = os.path.join(self.frequency_tif_path, "combined_frequency_floodmap.tif")
        with rasterio.open(combined_path, 'w', **meta) as dst:
            dst.write(combined_frequency, 1)
    
        print(f"Combined frequency map created at {combined_path}")


    def make_matrice_intensite(self):
        """
        Creates continuous polygons from valid pixels of the maximum depth GeoTIFF:
        - Level 1: Pixels with values > 0.3
        - Level 2: Pixels between 0 and 0.3
        
        No polygons are generated for nodata values.
        The result is saved into a new shapefile with fields:
        - "depth_zone": Description of the depth range.
        - "level": Integer representing the level (1 or 2).
        """
    
        # Define file paths
        depth_raster_path = os.path.join(self.max_depth_path, "combined_maximum_depth_map.tif")
        output_vector = os.path.join(self.max_depth_path, "split_depth_zones.shp")
    
        try:
            with rasterio.open(depth_raster_path) as src:
                raster_data = src.read(1)
                nodata = src.nodata
                transform = src.transform
                crs = src.crs
    
                # Create list for polygons
                features = []
    
                # Extract polygons for pixels > 0.3
                deep_mask = (raster_data > 0.3) & (raster_data != nodata)
                for geom, value in shapes(raster_data, mask=deep_mask, transform=transform):
                    poly = shape(geom)
                    features.append({
                        "properties": {
                            "depth_zone": "greater_than_0.3",
                            "level": 1
                        },
                        "geometry": poly
                    })
    
                # Extract polygons for pixels between 0 and 0.3
                shallow_mask = (raster_data <= 0.3) & (raster_data > 0) & (raster_data != nodata)
                for geom, value in shapes(raster_data, mask=shallow_mask, transform=transform):
                    poly = shape(geom)
                    features.append({
                        "properties": {
                            "depth_zone": "between_0_and_0.3",
                            "level": 2
                        },
                        "geometry": poly
                    })
    
                # Create a GeoDataFrame and save to a shapefile
                gdf = gpd.GeoDataFrame.from_features(features, crs=crs)
                gdf.to_file(output_vector, driver="ESRI Shapefile")
                print(f"Polygons saved to {output_vector}")
    
        except Exception as e:
            print(f"An error occurred: {e}")

    def process_objectives_de_protection(self):
        """
        Processes the maximum WSE raster by:
        1. Reclassifying values by rounding up to the first decimal place.
        2. Saving the reclassified raster.
        3. Extracting polygons from the reclassified raster.
        4. Merging small polygons (< 100 m²) into larger neighbors.
        5. Removing polygons corresponding to NoData values.
        """
    
        # Define file paths
        combined_path = os.path.join(self.max_wse_path, "combined_maximum_wse_map.tif")
        output_vector = os.path.join(self.obj_prot_path, "objectives_de_protection.shp")
        output_reclassified_raster = os.path.join(self.obj_prot_path, "objectives_de_protection.tif")
    
        # Open the input raster
        with rasterio.open(combined_path) as src:
            raster_data = src.read(1)
            profile = src.profile
            nodata = src.nodata if src.nodata is not None else -9999.0
    
            # Mask NoData values
            raster_data = np.where(raster_data == nodata, np.nan, raster_data)
    
            # Reclassify by rounding up to one decimal place
            reclassified_data = np.ceil(raster_data * 10) / 10
    
            # Replace NaN values with the original NoData value
            reclassified_data = np.where(np.isnan(reclassified_data), nodata, reclassified_data)
    
            # Ensure float32 compatibility
            reclassified_data = reclassified_data.astype(np.float32)
    
            # Update metadata and save the reclassified raster
            profile.update(dtype=rasterio.float32, count=1, nodata=nodata)
    
            with rasterio.open(output_reclassified_raster, 'w', **profile) as dst:
                dst.write(reclassified_data, 1)
    
            print("Generated protection objectives geotiff.")
    
            # Extract valid polygons
            mask = (reclassified_data != nodata)
            results = (
                {
                    "properties": {"wse_value": format(value, ".1f")},
                    "geometry": shape(geom),
                }
                for geom, value in shapes(reclassified_data, mask=mask, transform=src.transform)
                if value != nodata  # Ensure we exclude NoData polygons
            )
    
            # Create GeoDataFrame and filter out NoData values
            gdf = gpd.GeoDataFrame.from_features(results, crs=src.crs)
            gdf = gdf[gdf["wse_value"] != str(nodata)]  # Remove NoData polygons
    
        # Identify small polygons
        small_polygons = gdf[gdf.geometry.area < 100]
    
        # Process small polygons
        for idx, small_poly in small_polygons.iterrows():
            neighbors = gdf[gdf.geometry.touches(small_poly.geometry)]
    
            if not neighbors.empty:
                # Find the largest shared boundary
                best_match_idx = None
                max_shared_length = 0
    
                for neighbor_idx, neighbor in neighbors.iterrows():
                    shared_boundary = small_poly.geometry.intersection(neighbor.geometry).length
    
                    if shared_boundary > max_shared_length:
                        max_shared_length = shared_boundary
                        best_match_idx = neighbor_idx
    
                # Merge small polygon into the best-matching neighbor
                if best_match_idx is not None:
                    gdf.loc[best_match_idx, 'geometry'] = unary_union(
                        [gdf.loc[best_match_idx, 'geometry'], small_poly.geometry]
                    )
                    gdf.drop(idx, inplace=True)
    
        # Save the updated shapefile
        gdf.to_file(output_vector, driver="ESRI Shapefile")
        print("Generated protection objectives shapefile.")


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
    
        # Delete all folders matching folder_name_<number>
        for folder in os.listdir(parent_dir):
            if re.match(rf"{re.escape(folder_name)}_\d+$", folder):
                folder_to_remove = os.path.join(parent_dir, folder)
                if os.path.exists(folder_to_remove):
                    shutil.rmtree(folder_to_remove)
        
        # Delete all scripts matching script_base_name_<number>.py
        for script in os.listdir(script_dir):
            if re.match(rf"{re.escape(script_base_name)}_\d+\.py$", script):
                script_to_remove = os.path.join(script_dir, script)
                if os.path.exists(script_to_remove):
                    os.remove(script_to_remove)

        script_to_remove = os.path.join(script_dir, f"{script_base_name}.py")
        if os.path.exists(script_to_remove):
            os.remove(script_to_remove)

        
        print("Cleanup of parallel processing folders and scripts completed.")



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
        input_path = os.path.join(self.wse_profiles_path)

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
