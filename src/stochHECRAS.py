

import pandas as pd

import os
import glob
import shutil
import random
import time
import win32com.client
import numpy as np
import statistics
import rasterio
# import re
import csv
import h5py

from scipy.ndimage import label


class StochHECRAS:
    
        
    def __init__(self, stochICE):
        """
        Initializes the StochHECRAS class with stochastic input configurations.
    
        Args:
            stochICE (object): Object containing input configurations for the stochastic simulation.
        """
        self.stochICE = stochICE
        self.RC = None  # HEC-RAS Controller instance
    
        # Check if seed_parameters.csv exists in the project directory
        self.seed_file_path = os.path.join(self.stochICE.prjDir, "seed_parameters.csv")
        
        if os.path.exists(self.seed_file_path):
            print("\n'seed_parameters.csv' file found. Values of specified variables will be reused.")
        else:
            print("\nAttention: No 'seed_parameters.csv' file found. New stochastic variables will be generated.")
        
        self._load_seed_parameters()
    
    def _load_seed_parameters(self):
        if not os.path.exists(self.seed_file_path):
            self.parameters = []
            self.seed_data = None
            return
        
        with open(self.seed_file_path, "r") as file:
            lines = file.readlines()
        
        self.parameters = [line.strip() for line in lines if "," not in line and line.strip()]
        self.seed_data = pd.read_csv(self.seed_file_path, skiprows=len(self.parameters), header=0)
    
    def extract_or_randomize_data(self, index):
        if self.seed_data is not None and index < len(self.seed_data):
            row = self.seed_data.iloc[index]
        else:
            row = None
        
        mappings = {
            "Q": "flow",
            "phi": "phi",
            "porosity": "porosity",
            "Frontthick": "ice_thickness",
            'ice_cover_n':"ice_cover_n",
            "jam_loc_upstream": "location_upstream",
            "jam_loc_downstream": "location_downstream"
        }
        
        for param, attr in mappings.items():
            if param in self.parameters:
                setattr(self, attr, row[param] if row is not None else self._randomize_variable(attr))
            else:
                setattr(self, attr, self._randomize_variable(attr))
        
        if "jam_loc_upstream" in self.parameters and "jam_loc_downstream" in self.parameters:
            self.location = [self.location_downstream, self.location_upstream]
        else:
            self.location = random.choice(self.stochICE.jam_locations)
        
        self.apply_settings()
    
    def _randomize_variable(self, var_name):
        random_mappings = {
            "flow": random.choice(self.stochICE.Q),
            "phi": random.choice(self.stochICE.friction_angle),
            "porosity": random.choice(self.stochICE.porosity),
            "ice_thickness": random.choice(self.stochICE.Frontthick),
            "ice_cover_n": random.choice(self.stochICE.ice_cover_n),
            "ds_elev":random.choice(self.stochICE.ds_elev),
        }
        return random_mappings.get(var_name, None)
    
    def apply_settings(self):
        """
        Ensures that parameters are properly applied in the system.
        """
        self.set_flowrate()
        self.set_phi()
        self.set_porosity()
        self.thicknesses = [[f"{self.ice_thickness},{self.ice_thickness},{self.ice_thickness}"]]
        self.set_init_ice_cover_thickness()
        self.ice_cover_Manning = [[f"{self.ice_cover_n},{self.ice_cover_n},{self.ice_cover_n}"]]
        self.set_ice_cover_manning_values()
        self.set_ice_jam_location()

    def __getstate__(self):
        """
        Prepares the object state for pickling by removing unpicklable components.

        Removes the HEC-RAS controller object and other unpicklable attributes.

        Returns:
            dict: The object's state without unpicklable components.
        """
        state = self.__dict__.copy()
        state.pop("RC", None)  # Remove HEC-RAS Controller object
        state.pop("read_geofile", None)  # Remove other unpicklable attributes
        state.pop("out", None)
        return state

    def __setstate__(self, state):
        """
        Restores the object state after unpickling, reinitializing removed components.

        Args:
            state (dict): The state to restore.
        """
        self.__dict__.update(state)
        self.RC = win32com.client.Dispatch("RAS641.HECRASCONTROLLER")

    def preprocess_sims(self):
        """
        Prepares geometry nodes and initializes simulation parameters.

        This method opens the project file using the HEC-RAS controller and retrieves
        geometry nodes and other relevant data for the simulation.
        """
        self.RC = win32com.client.Dispatch("RAS641.HECRASCONTROLLER")
        self.RC.Project_Open(self.stochICE.ras_file)

        # Initialize geometry node variables
        self.NNod = None
        self.TabRS = None
        self.TabNTyp = None
        self.RiverID = 1
        self.ReachID = 1

        # Retrieve geometry nodes
        self.V1, self.v2, self.NNod, self.TabRS, self.TabNTyp = self.RC.Geometry_GetNodes(
            self.RiverID, self.ReachID, self.NNod, self.TabRS, self.TabNTyp
        )

        self.RC.QuitRAS()

        # Set output variable IDs
        self.WSE_id = 2  # Water Surface Elevation
        self.ice_thickness_id = 184  # Ice thickness ID
        self.MinChEle = 5  # Minimum channel elevation
        self.station_nbr_id = 161  # Station number ID

        # Store modified cross-section data
        self.xs_data_modified = self.stochICE.xs_data

    def new_launch_sims(self):
        """
        Launches stochastic simulations. For each simulation, it modifies parameters, writes new files, 
        and runs HEC-RAS to compute results for water surface elevation, ice thickness, and other variables.
        """
        print('\n-----------------------------------------------------------')
        print('                        Simulating')
        print('-----------------------------------------------------------\n')
    
        times = []
        self.result_profiles = {}
        self.input_parms = pd.DataFrame()
    
        self.sim_keys = []
        self.flow_rates = []
        self.init_ice_thicknesses = []
        self.phi_values = []
        self.porosity_values=[]
        self.ice_cover_n_values=[]
        self.jam_loc_upstream = []
        self.jam_loc_downstream = []
    
        # Read the initial geometry file content
        self.get_init_geofile_content()
    
        for j in range(self.stochICE.NSims):
            self.sim_key = f'sim_{j+1}'
            self.result_profiles[self.sim_key] = {}
    
            # Start the stopwatch for this simulation
            stopwatch = Stopwatch()
            stopwatch.start()
    
            # Randomize parameters and set up the simulation
            self.extract_or_randomize_data(j)

            # self.randomize_variables()
            self.write_new_geometry()
    
            # Store input variables
            self.sim_keys.append(self.sim_key)
            self.flow_rates.append(self.flow)
            self.init_ice_thicknesses.append(self.ice_thickness)
            self.phi_values.append(self.phi)
            self.porosity_values.append(self.porosity)
            self.ice_cover_n_values.append(self.ice_cover_n)
            self.jam_loc_upstream.append(self.location[1])
            self.jam_loc_downstream.append(self.location[0])
    
            # Save copies of geo and flow files
            self.store_geofile()
            self.store_flowfile()
    
            # Open the HEC-RAS project and run the simulation
            self.RC.Project_Open(self.stochICE.ras_file)
            self.NMsg, self.TabMsg, self.block = None, None, True
    
            print(f'Sim {j+1}/{self.stochICE.NSims} | Q = {self.flow:.2f}, Phi = {self.phi:.2f}, Porosity = {self.porosity:.2f}, Ice thickness = {self.ice_thickness:.2f}, Ice Mann = {self.ice_cover_n:.3f}')
            print(f"Ice Jam: Chainage {self.location[0]} - {self.location[1]}")
    
            # Run the current simulation plan
            self.v1, self.NMsg, self.TabMsg, self.v2 = self.RC.Compute_CurrentPlan(self.NMsg, self.TabMsg, self.block)
    
            # Collect simulation data
            wse_list, ice_thick_list, min_channel_ele_list, stations = [], [], [], []
    
            self.V1, self.v2, self.NNod, self.TabRS, self.TabNTyp = self.RC.Geometry_GetNodes(
                self.RiverID, self.ReachID, self.NNod, self.TabRS, self.TabNTyp
            )
    
            for i in range(self.NNod):
                if self.TabNTyp[i] == "":
                    # Water Surface Elevation
                    wse, *_ = self.RC.Output_NodeOutput(self.RiverID, self.ReachID, i + 1, 0, 1, self.WSE_id)
                    wse_list.append(wse)
    
                    # Ice Thickness
                    thick, *_ = self.RC.Output_NodeOutput(self.RiverID, self.ReachID, i + 1, 0, 1, self.ice_thickness_id)
                    ice_thick_list.append(thick)
    
                    # Minimum channel elevation
                    minChEle, *_ = self.RC.Output_NodeOutput(self.RiverID, self.ReachID, i + 1, 0, 1, self.MinChEle)
                    min_channel_ele_list.append(minChEle)
    
                    # Station number
                    stationNbr, *_ = self.RC.Output_NodeOutput(self.RiverID, self.ReachID, i + 1, 0, 1, self.station_nbr_id)
                    stations.append(stationNbr)
    
            # Store results for this simulation
            self.result_profiles[self.sim_key] = {
                'WSE': np.asarray(wse_list),
                'IceThick': np.asarray(ice_thick_list),
                'MinChEle': np.asarray(min_channel_ele_list),
                'TopIceMaxDepth': np.asarray(wse_list) - np.asarray(min_channel_ele_list),
                'station': np.asarray(stations),
            }

            #WSE map
            self.retain_largest_pixel_area(self.stochICE.wse_map_path)

            tif_filename = os.path.join(
                self.stochICE.wse_tifs_path,
                f"WSE_{self.flow}_{self.ice_thickness}_{self.ice_cover_n}_{self.phi}_{self.porosity}_{self.location[0]}_{self.location[1]}.tif",
            )
            shutil.copyfile(self.stochICE.wse_map_path, tif_filename)
    
            # Remove temporary files
            for _file in glob.glob(os.path.join(self.stochICE.wse_tifs_path, "*.vrt")):
                os.remove(_file)
    
    
            #Depth surface map
            self.retain_largest_pixel_area(self.stochICE.depth_map_path)

            tif_filename = os.path.join(
                self.stochICE.depth_tifs_path,
                f"Depth_{self.flow}_{self.ice_thickness}_{self.ice_cover_n}_{self.phi}_{self.porosity}_{self.location[0]}_{self.location[1]}.tif",
            )
            shutil.copyfile(self.stochICE.depth_map_path, tif_filename)
    
            # Remove temporary files
            for _file in glob.glob(os.path.join(self.stochICE.depth_tifs_path, "*.vrt")):
                os.remove(_file)

            self.extract_wse_profiles()
    
            # Stop the stopwatch and calculate simulation time
            sim_time = stopwatch.stop()
            times.append(sim_time)
            avg_time = statistics.mean(times)
            remaining_time = (self.stochICE.NSims - (j + 1)) * avg_time / 60
            print(f'Sim time: {sim_time:.2f}s, Avg time: {avg_time:.2f}s')
            print(f"Estimated time remaining: {remaining_time:.2f} mins\n")
    
        self.RC.QuitRAS()
    
        # Store input parameters across simulations
        self.input_parms = pd.DataFrame({
            #'sim_key': self.sim_keys,
            'Q': self.flow_rates,
            'phi': self.phi_values,
            'porosity': self.porosity_values,
            'Frontthick': self.init_ice_thicknesses,
            'ice_cover_n': self.ice_cover_n_values,
            'jam_loc_upstream': self.jam_loc_upstream,
            'jam_loc_downstream': self.jam_loc_downstream,
        })
        
        # Save to CSV
        csv_path = os.path.join(self.stochICE.prjDir, "sim_parameters.csv")
        self.input_parms.to_csv(csv_path, index=False)
        print(f"Simulation parameters saved to {csv_path}")          
        print(f"\n                  Simulations finished!\n")

    def retain_largest_pixel_area(self, raster_path):
        """
        Retains only the largest contiguous area of valid pixels in the specified raster file.
        Overwrites the original raster file with the result.
        
        This is useful to remove the 'puddles' not connected to the river left by HECRAS.
    
        Args:
            self: Instance of the class containing the function.
            raster_path (str): Path to the raster file to process.
        """

        # Load the raster
        with rasterio.open(raster_path) as src:
            raster_data = src.read(1)  # Read the first band
            meta = src.meta  # Metadata for the raster
            nodata_value = src.nodata  # Fetch the nodata value from the metadata
    
        # Identify valid pixels (non-nodata)
        valid_mask = raster_data != nodata_value
    
        # Label connected components of valid pixels
        labeled_array, num_features = label(valid_mask)
    
        # Calculate the size of each labeled region
        region_sizes = np.bincount(labeled_array.ravel())
    
        # Ignore the background region (label 0)
        region_sizes[0] = 0
    
        # Find the label of the largest region
        largest_region_label = region_sizes.argmax()
    
        # Create a mask for the largest region
        largest_region_mask = (labeled_array == largest_region_label)
    
        # Prepare the output data: Keep original values only in the largest region
        output_data = np.where(largest_region_mask, raster_data, nodata_value)
    
        # Update metadata for writing the output raster
        meta.update(dtype=rasterio.float32, compress='DEFLATE')
    
        # Overwrite the original file with the processed data
        with rasterio.open(raster_path, 'w', **meta) as dst:
            dst.write(output_data.astype(rasterio.float32), 1)
    
    def extract_wse_profiles(self):
        """
        Extracts water surface elevation (WSE) profiles for all unique reaches in the HDF5 file
        and saves them in the specified folder structure.
    
        Returns:
            None
        """
        plan_path = None
    
        # Locate the HDF5 file
        for root, _, files in os.walk(self.stochICE.prjDir):
            for file in files:
                if file.endswith('.p01.hdf'):
                    plan_path = os.path.join(root, file)
                    break
            if plan_path:
                break
    
        if not plan_path:
            raise FileNotFoundError("HDF5 plan file not found.")
    
        with h5py.File(plan_path, "r") as f:
            cross = f['Results']['Steady']['Output']['Geometry Info']['Cross Section Attributes'][()]
            wse = f['Results']['Steady']['Output']['Output Blocks']['Base Output']['Steady Profiles']['Cross Sections']['Water Surface'][()][0]
    
            reach_names = {entry[0].decode('utf-8') for entry in cross}
    
            for reach_name in reach_names:
                station_numbers, station_indices = [], []
    
                for i, entry in enumerate(cross):
                    if entry[0].decode('utf-8') == reach_name:
                        station_numbers.append(float(entry[2].decode('utf-8')))
                        station_indices.append(i)
    
                wse_values = [wse[idx] for idx in station_indices]
    
                wse_filename = os.path.join(
                    self.stochICE.wse_profiles_path,
                    f"{reach_name}_WSE_{self.flow}_{self.ice_thickness}_{self.ice_cover_n}_{self.phi}_{self.porosity}_{self.location[0]}_{self.location[1]}.csv"
                )
    
                with open(wse_filename, mode='w', newline='') as csv_file:
                    writer = csv.writer(csv_file)
                    writer.writerow(['Chainage (m)', 'wse (m)'])
                    for station, wse_value in zip(station_numbers, wse_values):
                        writer.writerow([station, wse_value])


    # def randomize_variables(self):
    #     """
    #     Randomizes key simulation parameters such as ice thickness, phi (friction angle), 
    #     and ice jam location based on specified ranges in `stochICE`.
    #     """
    #     self.xs_data_modified = self.stochICE.xs_data
    
    #     # Randomize flow rate
    #     self.flow = random.choice(self.stochICE.Q)
    #     self.set_flowrate()
    
    #     # Randomize phi (internal friction angle)
    #     self.phi = random.choice(self.stochICE.friction_angle)
    #     self.set_phi()
    
    #     # Randomize porosity
    #     self.porosity = random.choice(self.stochICE.porosity)
    #     self.set_porosity()
    
    #     # Randomize ice thickness
    #     self.ice_thickness = random.choice(self.stochICE.Frontthick)
    #     self.thicknesses = [[f"{self.ice_thickness},{self.ice_thickness},{self.ice_thickness}"]]
    #     self.set_init_ice_cover_thickness()
    
    #     # Randomize ice jam location
    #     self.location = random.choice(self.stochICE.jam_locations)
    #     self.set_ice_jam_location()
    
    def set_init_ice_cover_thickness(self):
        """Sets the initial ice thickness for each cross-section in the geometry."""
        for item in self.xs_data_modified:
            self.xs_data_modified[item]["Ice Thickness"]['val'] = self.thicknesses[0]

    def set_ice_cover_manning_values(self):
        """Sets the ice cover Manning values for each cross-section in the geometry."""
        for item in self.xs_data_modified:
            self.xs_data_modified[item]["Ice Mann"]['val'] = self.ice_cover_Manning[0]

    def set_phi(self):
        """Sets the ice friction angle for each cross-section in the geometry."""
        for item in self.xs_data_modified:
            self.xs_data_modified[item]['Ice Friction Angle']['val'] = self.phi
    
    def set_porosity(self):
        """Sets the porosity of the ice jam for each cross-section."""
        for item in self.xs_data_modified:
            self.xs_data_modified[item]['Ice Porosity']['val'] = self.porosity
    
    def set_ice_jam_location(self):
        """
        Updates the ice jam location for cross-sections where the reach matches 
        the ice jam reach and the chainage is within the specified range.
        """
        for item in self.xs_data_modified:
            if (self.xs_data_modified[item]['Reach'] == self.stochICE.ice_jam_reach and
                    self.location[0] <= self.xs_data_modified[item]['chainage'] <= self.location[1]):
                # Set ice jam location
                self.xs_data_modified[item]["Ice Is Channel"]['val'] = str(-1)
                self.xs_data_modified[item]["Ice Is OB"]['val'] = str(0)
            else:
                # Set default values for non-ice jam locations
                self.xs_data_modified[item]["Ice Is Channel"]['val'] = str(0)
                self.xs_data_modified[item]["Ice Is OB"]['val'] = str(0)
            
        self.remove_ice_upstream_of_jam()

    def remove_ice_upstream_of_jam(self):
        
        for item in self.xs_data_modified:
            if (self.xs_data_modified[item]['Reach'] == self.stochICE.ice_jam_reach and self.xs_data_modified[item]['chainage'] >= self.location[1]):
                self.xs_data_modified[item]["Ice Thickness"]['val'] = ["0,0,0"]
                
    
    def write_new_geometry(self):
        """
        Writes the updated geometry with modified ice cover parameters and ice jam locations 
        into the geometry file.
        """
        for key, item in self.xs_data_modified.items():
            # Update Ice Is Channel
            try:
                self.toWrite = f"Ice Is Channel={item['Ice Is Channel']['val']}\n"
                self.lineNmb = item['Ice Is Channel']['lnNum']
                self.replace_line()
            except TypeError:
                pass
    
            # Update Ice Is OB
            try:
                self.toWrite = f"Ice Is OB={item['Ice Is OB']['val']}\n"
                self.lineNmb = item['Ice Is OB']['lnNum']
                self.replace_line()
            except TypeError:
                pass
    
            # Update Ice Thickness
            try:
                self.toWrite = f"Ice Thickness={item['Ice Thickness']['val'][0]}\n"
                self.lineNmb = item['Ice Thickness']['lnNum']
                self.replace_line()
            except TypeError:
                pass

            # Update Ice Mann
            try:
                self.toWrite = f"Ice Mann={item['Ice Mann']['val'][0]}\n"
                self.lineNmb = item['Ice Mann']['lnNum']
                self.replace_line()
            except TypeError:
                pass

            # Update Ice Porosity
            try:
                
                porosity_value = item['Ice Porosity']['val']
                if isinstance(porosity_value, (list, tuple)):
                    porosity_value = porosity_value[0]
                self.toWrite = f"Ice Porosity={porosity_value}\n"
                # self.toWrite = f"Ice Porosity={item['Ice Porosity']['val'][0]}\n"
                self.lineNmb = item['Ice Porosity']['lnNum']
                self.replace_line()
            except TypeError:
                pass
    
            # Update Ice Friction Angle
            try:
                self.toWrite = f"Ice Friction Angle={item['Ice Friction Angle']['val']}\n"
                self.lineNmb = item['Ice Friction Angle']['lnNum']
                self.replace_line()
            except TypeError:
                pass
    
        # Finalize geometry file updates
        self.close_geofile()

    def get_init_geofile_content(self):
        """
        Reads the initial geometry file and stores its content.
    
        This method ensures that the original geometry file is loaded into memory
        for modifications during the stochastic simulations.
        """
        with open(self.stochICE.geo_file, 'r') as file:
            self.init_geofile_contents = file.readlines()
    
    def replace_line(self):
        """
        Replaces a specific line in the geometry file with updated content.
    
        This method modifies the in-memory representation of the geometry file
        and writes the changes back to the file.
        """
        self.new_geofile_contents = self.init_geofile_contents
        self.new_geofile_contents[self.lineNmb] = self.toWrite
    
        with open(self.stochICE.geo_file, 'w') as file:
            file.writelines(self.new_geofile_contents)
    
    def close_geofile(self):
        """
        Ensures the geometry file is properly closed.
    
        This method is mainly used to finalize and save the updated geometry file.
        """
        if hasattr(self, 'out') and self.out:
            self.out.close()

    def set_ds_elevation(self):
        """
        Sets downstream water level in the HEC-RAS flow file.

        """
        with open(self.stochICE.flow_file, 'r') as file:
            lines = file.readlines()
    
        updated_lines = []
        river_reaches = []
    
        # Identify the river reaches
        for i, line in enumerate(lines):
            if "Dn Known WS=" in line:
                river_reaches.append(i)
            updated_lines.append(line)
    
        # Update water level
        for idx, reach_line in enumerate(river_reaches):
            if idx == 0:
                updated_lines[reach_line + 1] = f"Dn Known WS={self.ds_elev}\n"
    
        # Save the updated flow file
        with open(self.stochICE.flow_file, 'w') as file:
            file.writelines(updated_lines)    



    def set_flowrate(self):
        """
        Updates the flow rate in the HEC-RAS flow file.
    
        This method supports scenarios with a single reach or an upstream reach
        divided into two sections. Flow rates are either applied directly or split
        evenly between upstream sections.
        """
        with open(self.stochICE.flow_file, 'r') as file:
            lines = file.readlines()
    
        updated_lines = []
        river_reaches = []
        current_reach_index = -1
    
        # Identify the river reaches
        for i, line in enumerate(lines):
            if "River Rch & RM" in line:
                current_reach_index += 1
                river_reaches.append(i)
            updated_lines.append(line)
    
        # Update flow rates
        for idx, reach_line in enumerate(river_reaches):
            if idx == 0:
                updated_lines[reach_line + 1] = f"     {self.flow}\n"
            else:
                half_value = self.flow / 2
                updated_lines[reach_line + 1] = f"     {half_value}\n"
    
        # Save the updated flow file
        with open(self.stochICE.flow_file, 'w') as file:
            file.writelines(updated_lines)
    
    def store_geofile(self):
        """
        Saves a copy of the modified geometry file with unique simulation parameters.
    
        This method ensures each simulation run has its own geometry file for
        future reference and debugging purposes.
        """
        geo_copy_path = os.path.join(
            self.stochICE.data_geo_files_path,
            f"{self.flow}_{self.ice_thickness}_{self.ice_cover_n}_{self.phi}_{self.location[0]}_{self.location[1]}.g01"
        )
        shutil.copyfile(self.stochICE.geo_file, geo_copy_path)
    
    def store_flowfile(self):
        """
        Saves a copy of the modified flow file with unique simulation parameters.
    
        Similar to `store_geofile`, this method creates a unique flow file
        for each simulation to track parameter variations and results.
        """
        flow_copy_path = os.path.join(
            self.stochICE.data_flow_files_path,
            f"{self.flow}_{self.ice_thickness}_{self.ice_cover_n}_{self.phi}_{self.location[0]}_{self.location[1]}.f01"
        )
        shutil.copyfile(self.stochICE.flow_file, flow_copy_path)

    def make_frequency_flood_map(self):
        """
        Generates an ensemble frequency flood map by aggregating flood maps from multiple simulations.
        The output is saved as a GeoTIFF file with pixel values representing the percentage of inundation.
    
        Args:
            None
    
        Returns:
            None
        """
        start_time = time.time()
        
        # Get list of flood maps
        flood_maps = os.listdir(self.stochICE.depth_tifs_path)
        aggregated_map = None
        
        # Process each flood map
        for flood_map in flood_maps:
            map_path = os.path.join(self.stochICE.depth_tifs_path, flood_map)
            with rasterio.open(map_path) as tiff:
                arr = tiff.read(1)  # Read the first band (32-bit float)
                arr = np.where(arr > 0.001, 1, 0)  # Binary: inundated (1) or not (0), matching original threshold
    
                if aggregated_map is None:
                    aggregated_map = np.zeros(arr.shape, dtype=np.float32)  # Initialize as float32
                    aggregated_map += arr
                    tiff_profile = tiff.profile  # Store profile for output
                else:
                    aggregated_map += arr
    
        # Calculate the frequency of inundation as a percentage
        frequency_map = (aggregated_map / self.stochICE.NSims) * 100  # Keep as float32
    
        # Define output path for the frequency flood map
        floodmap_output = os.path.join(self.stochICE.frequency_tif_path, f"frequency_floodmap_{self.stochICE.ID}.tif")
        
        # Update profile for float32 output
        tiff_profile.update(dtype='float32', nodata=-9999)  # Preserve precision with nodata value
    
        # Write the frequency map to a GeoTIFF file
        with rasterio.open(
            floodmap_output, 'w', 
            driver='GTiff',
            height=aggregated_map.shape[0], 
            width=aggregated_map.shape[1],
            count=1, 
            dtype='float32',  # Use float32 instead of uint8
            crs=tiff_profile['crs'], 
            transform=tiff_profile['transform'],
            compress='deflate',  # Match original compression
            nodata=-9999
        ) as result:
            result.write(frequency_map, 1)
        
        elapsed_time = time.time() - start_time
        print(f"Frequency flood map created successfully in {elapsed_time:.2f} seconds.")

    
    def make_maximum_depth_map(self):
        """
        Generates a maximum depth map by comparing water depths from multiple GeoTIFF files.
        The output is saved as a GeoTIFF file with pixel values representing the maximum depth.
    
        Args:
            None
    
        Returns:
            None
        """
        start_time = time.time()  # Start timing the process
    
        depth_maps = os.listdir(self.stochICE.depth_tifs_path)
        max_depth = None
    
        for depth_map in depth_maps:
            map_path = os.path.join(self.stochICE.depth_tifs_path, depth_map)
            with rasterio.open(map_path) as tiff:
                arr = tiff.read(1)  # Read the first band
                
                if max_depth is None:
                    max_depth = arr
                else:
                    max_depth = np.maximum(max_depth, arr)
    
        # Define output path for the maximum depth map
        max_depth_output = os.path.join(self.stochICE.max_depth_tif_path, f"maximum_depth_map_{self.stochICE.ID}.tif")
    
        # Write the maximum depth map to a GeoTIFF file
        with rasterio.open(
            max_depth_output, 'w', driver='GTiff',
            height=max_depth.shape[0], width=max_depth.shape[1],
            count=1, dtype=max_depth.dtype, crs=tiff.crs, transform=tiff.transform,
            compress='deflate'  # Use Deflate compression for smaller file size
        ) as result:
            result.write(max_depth, 1)
    
        elapsed_time = time.time() - start_time  # Calculate elapsed time
        print(f"Maximum depth map created successfully in {elapsed_time:.2f} seconds.")


    def make_maximum_wse_map(self):
        """
        Generates a maximum wse map by comparing water levels of the individual GeoTIFF files.
        The output is saved as a GeoTIFF file with pixel values representing the maximum wse.
    
        Args:
            None
    
        Returns:
            None
        """
        start_time = time.time()  # Start timing the process
    
        wse_maps = os.listdir(self.stochICE.wse_tifs_path)
        max_wse = None
    
        for wse_map in wse_maps:
            map_path = os.path.join(self.stochICE.wse_tifs_path, wse_map)
            with rasterio.open(map_path) as tiff:
                arr = tiff.read(1)  # Read the first band
                
                if max_wse is None:
                    max_wse = arr
                else:
                    max_wse = np.maximum(max_wse, arr)
    
        # Define output path for the maximum depth map
        max_wse_output = os.path.join(self.stochICE.max_wse_tif_path, f"maximum_wse_map_{self.stochICE.ID}.tif")
    
        # Write the maximum depth map to a GeoTIFF file
        with rasterio.open(
            max_wse_output, 'w', driver='GTiff',
            height=max_wse.shape[0], width=max_wse.shape[1],
            count=1, dtype=max_wse.dtype, crs=tiff.crs, transform=tiff.transform,
            compress='deflate'  # Use Deflate compression for smaller file size
        ) as result:
            result.write(max_wse, 1)
    
        elapsed_time = time.time() - start_time  # Calculate elapsed time
        print(f"Maximum wse map created successfully in {elapsed_time:.2f} seconds.")



class Stopwatch:
    """
    A simple stopwatch class to measure elapsed time for simulations.
    """
    def __init__(self):
        self.start_time = None

    def start(self):
        """Starts the stopwatch."""
        self.start_time = time.time()

    def stop(self):
        """Stops the stopwatch and returns the elapsed time."""
        if self.start_time is None:
            raise ValueError("Stopwatch has not been started.")
        elapsed_time = time.time() - self.start_time
        self.start_time = None
        return elapsed_time



# class StochHECRAS:
#     """
#     This class manages stochastic simulations using the HEC-RAS controller. 
#     It handles the preparation, execution, and post-processing of multiple simulations 
#     by modifying key input parameters like flow rates and ice thicknesses.
#     """
    
#     def __init__(self, stochICE):
#         """
#         Initialize the StochHECRAS class with stochastic input configuration.

#         Args:
#             stochICE: An object containing input configurations for the stochastic simulation.
#         """
#         self.stochICE = stochICE

#         # Optionally clear previous results
#         # if self.stochICE.clr:
#         #     self.clear_results()

#     def __getstate__(self):
#         """
#         Prepares the object state for pickling by removing certain components 
#         that should not be pickled (e.g., the HECRAS controller object).
#         """
#         state = self.__dict__.copy()
#         del state["RC"]
#         del state["read_geofile"]
#         del state["out"]
#         return state

#     def __setstate__(self, state):
#         """
#         Restores the object state after unpickling, reinitializing components 
#         that were removed during pickling.
#         """
#         self.__dict__.update(state)
#         self.RC = win32com.client.Dispatch("RAS641.HECRASCONTROLLER")

#     def preprocess_sims(self):
#         """
#         Prepares the geometry nodes and initializes simulation parameters.
#         It also opens the project file using the HEC-RAS controller.
#         """
#         self.RC = win32com.client.Dispatch("RAS641.HECRASCONTROLLER")
#         self.RC.Project_Open(self.stochICE.ras_file)
#         self.NNod, self.TabRS, self.TabNTyp = None, None, None
#         self.RiverID = 1
#         self.ReachID = 1
#         # Get geometry nodes
#         self.V1, self.v2, self.NNod, self.TabRS, self.TabNTyp = self.RC.Geometry_GetNodes(
#             self.RiverID, self.ReachID, self.NNod, self.TabRS, self.TabNTyp)
#         self.RC.QuitRAS()

#         # Set output variable IDs
#         self.WSE_id = 2              # Water Surface Elevation
#         self.ice_thickness_id = 184   # Ice thickness ID
#         self.MinChEle = 5             # Minimum channel elevation
#         self.station_nbr_id = 161     # Station number ID

#         self.xs_data_modified = self.stochICE.xs_data

    # def new_launch_sims(self):
    #     """
    #     Launches stochastic simulations. For each simulation, it modifies parameters, writes new files, 
    #     and runs HEC-RAS to compute results for water surface elevation, ice thickness, and other variables.
    #     """
    #     print('-----------------------------------------------------------')
    #     print('-----------------------Simulating--------------------------')
    #     print('-----------------------------------------------------------\n')

    #     times = []
    #     self.result_profiles = {}
    #     self.input_parms = pd.DataFrame()

    #     self.sim_keys = []
    #     self.flow_rates = []
    #     self.init_ice_thicknesses = []
    #     self.phi_values = []
    #     self.jam_loc_upstream = []
    #     self.jam_loc_downstream = []

    #     # Read the initial geometry file content
    #     self.get_init_geofile_content()

    #     for j in range(self.stochICE.NSims):
    #         self.sim_key = f'sim_{j+1}'
    #         self.result_profiles[self.sim_key] = {}

    #         # Start the stopwatch for this simulation
    #         stopwatch = Stopwatch()
    #         stopwatch.start()

    #         # Randomize parameters and set up the simulation
    #         self.randomize_variables()
            
            
    #         self.write_new_geometry()
            
    #         # Store input variables
    #         self.sim_keys.append(self.sim_key)
    #         self.flow_rates.append(self.flow)
    #         self.init_ice_thicknesses.append(self.ice_thickness)
    #         self.phi_values.append(self.phi)
    #         self.jam_loc_upstream.append(self.location[0])
    #         self.jam_loc_downstream.append(self.location[1])

    #         # Save copies of geo and flow files
    #         self.store_geofile()
    #         self.store_flowfile()

    #         # Open the HEC-RAS project and run the simulation
    #         self.RC.Project_Open(self.stochICE.ras_file)
    #         self.NMsg, self.TabMsg, self.block = None, None, True

    #         print(f'Sim {j+1} of {self.stochICE.NSims}.')
    #         print(f"Running Q = {self.flow:.2f}, Phi = {self.phi:.2f}, Ice thickness = {self.ice_thickness:.2f}.")
    #         print(f"Ice Jam location between chainage {self.location[0]} and {self.location[1]}.")

    #         # Run the current simulation plan
    #         self.v1, self.NMsg, self.TabMsg, self.v2 = self.RC.Compute_CurrentPlan(self.NMsg, self.TabMsg, self.block)

    #         # Collect simulation data
    #         wse_list, ice_thick_list, min_channel_ele_list, stations = [], [], [], []

    #         # Get nodes' output data
    #         self.V1, self.v2, self.NNod, self.TabRS, self.TabNTyp = self.RC.Geometry_GetNodes(
    #             self.RiverID, self.ReachID, self.NNod, self.TabRS, self.TabNTyp)

    #         for i in range(self.NNod):
    #             if self.TabNTyp[i] == "":
    #                 # Water Surface Elevation
    #                 wse, *_ = self.RC.Output_NodeOutput(self.RiverID, self.ReachID, i+1, 0, 1, self.WSE_id)
    #                 wse_list.append(wse)

    #                 # Ice Thickness
    #                 thick, *_ = self.RC.Output_NodeOutput(self.RiverID, self.ReachID, i+1, 0, 1, self.ice_thickness_id)
    #                 ice_thick_list.append(thick)

    #                 # Minimum channel elevation
    #                 minChEle, *_ = self.RC.Output_NodeOutput(self.RiverID, self.ReachID, i+1, 0, 1, self.MinChEle)
    #                 min_channel_ele_list.append(minChEle)

    #                 # Station number
    #                 stationNbr, *_ = self.RC.Output_NodeOutput(self.RiverID, self.ReachID, i+1, 0, 1, self.station_nbr_id)
    #                 stations.append(stationNbr)

    #         # Store results for this simulation
    #         self.result_profiles[self.sim_key] = {
    #             'WSE': np.asarray(wse_list),
    #             'IceThick': np.asarray(ice_thick_list),
    #             'MinChEle': np.asarray(min_channel_ele_list),
    #             'TopIceMaxDepth': np.asarray(wse_list) - np.asarray(min_channel_ele_list),
    #             'station': np.asarray(stations),
    #         }

    #         # Save the 2D flood map
    #         # tif_filename = os.path.join(
    #         #     self.stochICE.prjDir, "MonteCarlo", "SimulationTifs",
    #         #     f"WSE_{self.flow}_{self.ice_thickness}_{self.phi}_{self.porosity}_{self.location[0]}_{self.location[1]}.tif"
    #         # )
            
    #         tif_filename = os.path.join(self.stochICE.depths_path,
    #             f"WSE_{self.flow}_{self.ice_thickness}_{self.phi}_{self.porosity}_{self.location[0]}_{self.location[1]}.tif"
    #         )            
            
    #         shutil.copyfile(self.stochICE.wse_map_path, tif_filename)

    #         # Remove temporary files
    #         # for _file in glob.glob(os.path.join(self.stochICE.prjDir, "MonteCarlo", "SimulationTifs", "*.vrt")):
    #         #     os.remove(_file)

    #         # Remove temporary files
    #         for _file in glob.glob(os.path.join(self.stochICE.depths_path, "*.vrt")):
    #             os.remove(_file)


    #         self.extract_wse_profiles()
    #         # Stop the stopwatch and calculate simulation time
    #         sim_time = stopwatch.stop()
    #         times.append(sim_time)
    #         avg_time = statistics.mean(times)
    #         print(f'Sim time {sim_time:.2f} s, average time {avg_time:.2f} s.')
    #         remaining_time = (self.stochICE.NSims - (j+1)) * avg_time / 60
    #         print(f"Approximately {remaining_time:.2f} mins remaining in batch.\n")

    #     self.RC.QuitRAS()

    #     # Store input parameters across simulations
    #     self.input_parms = pd.DataFrame({
    #         'sim_key': self.sim_keys,
    #         'Q': self.flow_rates,
    #         'phi': self.phi_values,
    #         'init_ice_thicknesses': self.init_ice_thicknesses,
    #         'jam_loc_upstream': self.jam_loc_upstream,
    #         'jam_loc_downstream': self.jam_loc_downstream,
    #     })

    # def extract_wse_profiles(self):
    #     """
    #     Detects all unique reaches in the file, and for each reach, extracts station numbers
    #     and water surface elevations, saving them to a respective CSV file in a specified folder structure.
        
    #     Parameters:
    #     filename (str): Path to the HDF5 file.
    #     output_folder (str): Path to the main output folder where subfolders will be created.
        
    #     Returns:
    #     None
    #     """
        
    #     for root, dirs, files in os.walk(self.stochICE.prjDir):
    #         for file in files:
    #             if file.endswith('.p01.hdf'):
    #                 plan_path = os.path.join(root, file)
                    
    #                 # return os.path.join(root, file)
    
    #     # r'C:\Users\dugj2403\Desktop\StochICE_refine_project\Chateauguay_ArthurLaberge\Chateauguay_JD.p01.hdf'
       
        
    #     # self.wse_path = os.path.join(self.stochICE.prjDir, "MonteCarlo", "WSE")
    #     self.wse_path = os.path.join(self.stochICE.wse_path)


    #     # if not os.path.exists(self.wse_path):
    #     #     os.makedirs(self.wse_path)
        
    #     with h5py.File(plan_path, "r") as f:
    #         # Access the 'Cross Section Attributes' and 'Water Surface' datasets
    #         cross = f['Results']['Steady']['Output']['Geometry Info']['Cross Section Attributes'][()]
    #         wse = f['Results']['Steady']['Output']['Output Blocks']['Base Output']['Steady Profiles']['Cross Sections']['Water Surface'][()][0]
            
    #         # Detect all unique reach names
    #         reach_names = set(entry[0].decode('utf-8') for entry in cross)
            
    #         # Process each reach individually
    #         for reach_name in reach_names:
    #             station_numbers = []
    #             station_indices = []
                
    #             # Collect station numbers and indices for the current reach
    #             for i, entry in enumerate(cross):
    #                 if entry[0].decode('utf-8') == reach_name:
    #                     station_numbers.append(int(entry[2].decode('utf-8')))
    #                     station_indices.append(i)
                
    #             # Extract water surface elevations for the current reach
    #             wse_values = [wse[idx] for idx in station_indices]
                
    #             # Define the subfolder and CSV file path for the current reach
    #             #reach_folder = os.path.join(self.wse_path, reach_name)
    #             #os.makedirs(reach_folder, exist_ok=True)
    #             wse_filename = os.path.join(self.wse_path,
    #                 f"{reach_name}_WSE_{self.flow}_{self.ice_thickness}_{self.phi}_{self.porosity}_{self.location[0]}_{self.location[1]}.csv"
    #             )
                
    #             # output_csv = os.path.join(reach_folder, f"WSE_{self.flow}_{self.ice_thickness}_{self.phi}_{self.porosity}_{self.location[0]}_{self.location[1]}_{reach_name}.csv")
                
    #             # Write the station numbers and WSE values to the CSV file
    #             with open(wse_filename, mode='w', newline='') as csv_file:
    #                 writer = csv.writer(csv_file)
    #                 writer.writerow(['Chainage (m)', 'wse (m)'])  # Header row
    #                 for station, wse_value in zip(station_numbers, wse_values):
    #                     writer.writerow([station, wse_value])




    # def randomize_variables(self):
    #     """
    #     Randomizes key simulation parameters such as ice thickness, phi (friction angle), 
    #     and ice jam location based on specified ranges in `stochICE`.
    #     """

    #     self.xs_data_modified = self.stochICE.xs_data
        
    #     #Flow rate
    #     self.flow=random.sample(self.stochICE.Q, 1)[0]
    #     self.set_flowrate()
       
    #     #phi (internal friction angle)
    #     self.phi=random.sample(self.stochICE.friction_angle, 1)[0]
    #     self.set_phi()
        
    #     #porosity 
    #     self.porosity=random.sample(self.stochICE.porosity, 1)[0]
    #     self.set_porosity()
        
    #     # thicknesses
    #     self.ice_thickness=random.sample(self.stochICE.Frontthick, 1)[0]
    #     self.thicknesses = [[f"{self.ice_thickness},{self.ice_thickness},{self.ice_thickness}"]]
    #     self.set_init_ice_cover_thickness()
        
    #     # Randomly select ice jam location from the available locations
    #     self.location = random.choice(self.stochICE.jam_locations)
    #     self.set_ice_jam_location()

    # def set_init_ice_cover_thickness(self):
    #     """Sets the initial ice thickness for each cross-section in the geometry."""
    #     for item in self.xs_data_modified:
    #         self.xs_data_modified[item]["Ice Thickness"]['val'] = self.thicknesses[0]

    # def set_phi(self):
    #     """Sets the ice friction angle for each cross-section in the geometry."""
    #     for item in self.xs_data_modified:
    #         self.xs_data_modified[item]['Ice Friction Angle']['val'] = self.phi

    # def set_porosity(self):
    #     """Sets the porosity of the ice jam."""
    #     for item in self.xs_data_modified:
    #         self.xs_data_modified[item]['Ice Porosity']['val'] = self.porosity

    # def set_ice_jam_location(self):
    #     """
    #     Updates the ice jam location for cross-sections where the reach matches 
    #     the ice jam reach and the chainage is within the specified range.
    #     """
    #     for item in self.xs_data_modified:
    #         if (self.xs_data_modified[item]['Reach'] == self.stochICE.ice_jam_reach and
    #                 self.xs_data_modified[item]['chainage'] >= self.location[0] and
    #                 self.xs_data_modified[item]['chainage'] <= self.location[1]):
    #             # Set ice jam location
    #             self.xs_data_modified[item]["Ice Is Channel"]['val'] = str(-1)
    #             self.xs_data_modified[item]["Ice Is OB"]['val'] = str(0)
    #         else:
    #             # Set default values for non-ice jam locations
    #             self.xs_data_modified[item]["Ice Is Channel"]['val'] = str(0)
    #             self.xs_data_modified[item]["Ice Is OB"]['val'] = str(0)

    # def write_new_geometry(self):
    #     """
    #     Writes the new geometry with modified ice cover and ice jam locations 
    #     into the geometry file.
    #     """
    #     for key, item in self.xs_data_modified.items():
    #         try:
    #             self.toWrite = f"Ice Is Channel={item['Ice Is Channel']['val']}\n"
    #             self.lineNmb = item['Ice Is Channel']['lnNum']
    #             self.replace_line()
    #         except TypeError:
    #             pass

    #         try:
    #             self.toWrite = f"Ice Is OB={item['Ice Is OB']['val']}\n"
    #             self.lineNmb = item['Ice Is OB']['lnNum']
    #             self.replace_line()
    #         except TypeError:
    #             pass

    #         try:
    #             self.toWrite = f"Ice Thickness={item['Ice Thickness']['val'][0]}\n"
    #             self.lineNmb = item['Ice Thickness']['lnNum']
    #             self.replace_line()
    #         except TypeError:
    #             pass

    #         try:
    #             self.toWrite = f"Ice Porosity={item['Ice Porosity']['val'][0]}\n"
    #             self.lineNmb = item['Ice Porosity']['lnNum']
    #             self.replace_line()
    #         except TypeError:
    #             pass

    #         try:
    #             self.toWrite = f"Ice Friction Angle={item['Ice Friction Angle']['val']}\n"
    #             self.lineNmb = item['Ice Friction Angle']['lnNum']
    #             self.replace_line()
    #         except TypeError:
    #             pass

    #     self.close_geofile()

    # def get_init_geofile_content(self):
    #     """Reads the initial geometry file and stores its content."""
    #     self.read_geofile = open(self.stochICE.geo_file, 'r')
    #     self.init_geofile_contents = self.read_geofile.readlines()
    #     self.read_geofile.close()

    # def replace_line(self):
    #     """Replaces a line in the geometry file with modified content."""
    #     self.new_geofile_contents = self.init_geofile_contents
    #     self.new_geofile_contents[self.lineNmb] = self.toWrite

    #     with open(self.stochICE.geo_file, 'w') as self.out:
    #         self.out.writelines(self.new_geofile_contents)

    # def close_geofile(self):
    #     """Closes the geometry file."""
    #     self.out.close()

    # def set_flowrate(self):
    #     """
    #     Updates the flow rate in the HEC-RAS flow file. 
    #     This implementation assumes a single reach or a divided upstream reach.
    #     """
    #     with open(self.stochICE.flow_file, 'r') as file:
    #         lines = file.readlines()

    #     updated_lines = []
    #     river_reaches = []
    #     current_reach_index = -1

    #     for i, line in enumerate(lines):
    #         if "River Rch & RM" in line:
    #             current_reach_index += 1
    #             river_reaches.append(i)
    #         updated_lines.append(line)

    #     for idx, reach_line in enumerate(river_reaches):
    #         if idx == 0:
    #             updated_lines[reach_line + 1] = f"     {self.flow}\n"
    #         else:
    #             half_value = self.flow / 2
    #             updated_lines[reach_line + 1] = f"     {half_value}\n"

    #     with open(self.stochICE.flow_file, 'w') as file:
    #         file.writelines(updated_lines)

    # def store_geofile(self):
    #     """Saves a copy of the modified geometry file with unique simulation parameters."""
    #     geo_copy_path = os.path.join(
    #         self.stochICE.data_geo_files_path,
    #         f"{self.flow}_{self.ice_thickness}_{self.phi}_{self.location[0]}_{self.location[1]}.g01")
    #     shutil.copyfile(self.stochICE.geo_file, geo_copy_path)

    # def store_flowfile(self):
    #     """Saves a copy of the modified flow file with unique simulation parameters."""
    #     flow_copy_path = os.path.join(
    #         self.stochICE.data_flow_files_path,
    #         f"{self.flow}_{self.ice_thickness}_{self.phi}_{self.location[0]}_{self.location[1]}.f01")
    #     shutil.copyfile(self.stochICE.flow_file, flow_copy_path)

    # def clear_results(self):
    #     """
    #     Deletes previous simulation results from the MonteCarlo folder, including *.tifs, 
    #     *.g0, and *.f0 files.
    #     """
    #     print('Deleting *.tifs, *.g0 and *.f0 records in MonteCarlo folder!\nFlag "clrRes=False" to suppress.\n')
    #     files = glob.glob(self.stochICE.tif_path + '\\*')
    #     for f in files:
    #         os.remove(f)

    #     files = glob.glob(self.stochICE.geo_files_path + '\\*')
    #     for f in files:
    #         os.remove(f)

    #     files = glob.glob(self.stochICE.flow_files_path + '\\*')
    #     for f in files:
    #         os.remove(f)

    # def compress_results(self):
    #     """Compresses the simulation results into a zip file."""
    #     shutil.make_archive(
    #         self.stochICE.prjDir + '\\MC_results_batch_%s' % self.stochICE.ID,
    #         'zip',
    #         self.stochICE.MC_path)

    # def make_frequency_flood_map(self):
    #     """
    #     Generates an ensemble flood map by aggregating the flood maps generated from 
    #     multiple simulations and saves it as a TIFF file.
    #     """
    #     listdir = os.listdir(self.stochICE.depths_path)

    #     for count, m in enumerate(listdir):
    #         data_name = os.path.join(self.stochICE.depths_path, m)
    #         tiff = rasterio.open(data_name)
    #         arr = tiff.read()

    #         arr[arr > 0] = 1
    #         arr[arr < 0] = 0

    #         if count == 0:
    #             stoch = arr
    #         else:
    #             stoch = stoch + arr

    #     self.stoch = (stoch / self.stochICE.NSims) * 100
    #     # self.floodmap_path = os.path.join(self.stochICE.frequency_path)

    #     # if not os.path.exists(self.floodmap_path):
    #     #     os.makedirs(self.floodmap_path)

    #     self.stoch = self.stoch.astype(int)

    #     # try:
    #     #     os.remove(os.path.join(self.stochICE.frequency_path, "frequency_floodmap.tif"))
    #     # except FileNotFoundError:
    #     #     pass

    #     floodmap_output = os.path.join(self.stochICE.frequency_path, f"frequency_floodmap_{self.stochICE.ID}.tif")

    #     result = rasterio.open(
    #         floodmap_output, 'w', driver='GTiff',
    #         height=tiff.shape[0], width=tiff.shape[1],
    #         count=1, dtype=stoch.dtype, crs=tiff.crs, transform=tiff.transform,
    #         compress='deflate')

    #     result.write(stoch[0, :, :], 1)
    #     result.close()
    #     print("Frequency map created.")



    # def make_maximum_depth_map(self):
    #     """
    #     Generates a maximum depth map by comparing the water depths 
    #     from multiple GeoTIFF files and saves it as a TIFF file.
    #     """
    
    #     # Start the timer
    #     start_time = time.time()
    
    #     listdir = os.listdir(self.stochICE.depths_path)
    
    #     # Initialize an array to store the maximum depths
    #     max_depth = None
    
    #     for count, m in enumerate(listdir):
    #         data_name = os.path.join(self.stochICE.depths_path, m)
    #         with rasterio.open(data_name) as tiff:
    #             arr = tiff.read(1)  # Read the first band
                
    #             if max_depth is None:
    #                 # Initialize max_depth with the same shape as arr
    #                 max_depth = arr
    #             else:
    #                 # Update max_depth array with the maximum values at each pixel
    #                 max_depth = np.maximum(max_depth, arr)
    
    #     # Define the path for the output maximum depth map
    #     # max_depth_map_path = os.path.join(self.stochICE.prjDir, "FloodMaps")
    #     # if not os.path.exists(max_depth_map_path):
    #     #     os.makedirs(max_depth_map_path)
    
    #     max_depth_output = os.path.join(self.stochICE.max_depth_path, f"maximum_depth_map_{self.stochICE.ID}.tif")
    
    #     # Save the maximum depth map as a GeoTIFF with Deflate compression
    #     with rasterio.open(
    #         max_depth_output, 'w', driver='GTiff',
    #         height=max_depth.shape[0], width=max_depth.shape[1],
    #         count=1, dtype=max_depth.dtype, crs=tiff.crs, transform=tiff.transform,
    #         compress='deflate'  # Use Deflate compression
    #     ) as result:
    #         result.write(max_depth, 1)
        
    #     # End the timer and calculate the elapsed time
    #     end_time = time.time()
    #     elapsed_time = end_time - start_time
    
    #     # Print the elapsed time
    #     print(f"Time required to calculate the maximum depth map: {elapsed_time:.2f} seconds.\n")


    # # def calculate_flood_areas(self):
    #     """
    #     Calculates the area of all pixels containing data in each GeoTIFF file, 
    #     extracts metadata from the file name, and saves the results to a CSV file.
    #     """
    
    #     # Start the timer
    #     start_time = time.time()
    
    #     listdir = os.listdir(self.stochICE.tif_path)
    #     results = []
    
    #     for m in listdir:
    #         data_name = os.path.join(self.stochICE.tif_path, m)
    
    #         # Extract data from the filename using regular expressions
    #         # Filename structure assumed to be: <Q>_<thick>_<phi>_<porosity>_<foot>_<head>.tif
    #         match = re.match(r"WSE_(\d+\.\d+)_(\d+\.\d+)_(\d+\.\d+)_(\d+\.\d+)_(\d+)_(\d+)", m)
    #         if match:
    #             Q, thick, phi, porosity, foot, head = match.groups()
    #         else:
    #             print(f"Filename '{m}' does not match the expected format.")
    #             continue
    
    #         with rasterio.open(data_name) as tiff:
    #             arr = tiff.read(1)  # Read the first band
    
    #             # Calculate the pixel area based on the spatial resolution
    #             pixel_area = abs(tiff.transform[0] * tiff.transform[4])  # Width * Height in projected units
    #             data_pixels = np.count_nonzero(arr != -9999)
    #             # data_pixels = np.count_nonzero(~np.isnan(arr))  # Count pixels with data (non-NaN)
    #             area_with_data = data_pixels * pixel_area  # Total area for pixels with data
    
    #             # Append result to list
    #             results.append({
    #                 "Q": Q,
    #                 "thick": thick,
    #                 "phi": phi,
    #                 "porosity": porosity,
    #                 "foot": foot,
    #                 "head": head,
    #                 "area_with_data": area_with_data
    #             })
    
    #     # Define the path for the output CSV file
    #     csv_output_path = os.path.join(self.stochICE.prjDir, "FloodMaps", "area_per_map.csv")
    
    #     # Ensure the directory exists
    #     os.makedirs(os.path.dirname(csv_output_path), exist_ok=True)
    
    #     # Write results to a CSV file
    #     with open(csv_output_path, mode='w', newline='') as csv_file:
    #         fieldnames = ["Q", "thick", "phi", "porosity", "foot", "head", "area_with_data"]
    #         writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
            
    #         writer.writeheader()
    #         for result in results:
    #             writer.writerow(result)
    
    #     # End the timer and calculate the elapsed time
    #     end_time = time.time()
    #     elapsed_time = end_time - start_time
    
    #     # Print the elapsed time
    #     print(f"Time required to calculate areas: {elapsed_time:.2f} seconds\n")


    # def perform_multivariate_regression(self):
    #     # Define the path to the file
    #     file_path = os.path.join(self.stochICE.prjDir, "FloodMaps", "area_per_map.csv")
        
    #     # Load the data
    #     data = pd.read_csv(file_path)
    
    #     # Define the dependent variable and independent variables
    #     dependent_var = 'area_with_data'
    #     independent_vars = data.columns.drop(dependent_var)
    
    #     # Separate the data into X (independent variables) and y (dependent variable)
    #     X = data[independent_vars]
    #     y = data[dependent_var]
    
    #     # Split the data into training and testing sets
    #     X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    
    #     # Initialize and fit the linear regression model
    #     model = LinearRegression()
    #     model.fit(X_train, y_train)
    
    #     # Predict on the test set
    #     y_pred = model.predict(X_test)
    
    #     # Calculate metrics
    #     mse = mean_squared_error(y_test, y_pred)
    #     r2 = r2_score(y_test, y_pred)
    
    #     # Display the coefficients
    #     coefficients = pd.Series(model.coef_, index=independent_vars)
    
    #     print("Regression Coefficients:")
    #     print(coefficients)
    #     print("\nModel Performance:")
    #     print(f"Mean Squared Error (MSE): {mse}")
    #     print(f"R-squared (R²): {r2}")
    
