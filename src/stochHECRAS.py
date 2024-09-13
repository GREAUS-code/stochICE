# -*- coding: utf-8 -*-
"""
Class of functions to run HECRAS stochastically
"""

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


class StochHECRAS:
    """
    This class manages stochastic simulations using the HEC-RAS controller. 
    It handles the preparation, execution, and post-processing of multiple simulations 
    by modifying key input parameters like flow rates and ice thicknesses.
    """
    
    def __init__(self, stochICE):
        """
        Initialize the StochHECRAS class with stochastic input configuration.

        Args:
            stochICE: An object containing input configurations for the stochastic simulation.
        """
        self.stochICE = stochICE

        # Optionally clear previous results
        if self.stochICE.clr:
            self.clear_results()

    def __getstate__(self):
        """
        Prepares the object state for pickling by removing certain components 
        that should not be pickled (e.g., the HECRAS controller object).
        """
        state = self.__dict__.copy()
        del state["RC"]
        del state["read_geofile"]
        del state["out"]
        return state

    def __setstate__(self, state):
        """
        Restores the object state after unpickling, reinitializing components 
        that were removed during pickling.
        """
        self.__dict__.update(state)
        self.RC = win32com.client.Dispatch("RAS641.HECRASCONTROLLER")

    def preprocess_sims(self):
        """
        Prepares the geometry nodes and initializes simulation parameters.
        It also opens the project file using the HEC-RAS controller.
        """
        self.RC = win32com.client.Dispatch("RAS641.HECRASCONTROLLER")
        self.RC.Project_Open(self.stochICE.ras_file)
        self.NNod, self.TabRS, self.TabNTyp = None, None, None
        self.RiverID = 1
        self.ReachID = 1
        # Get geometry nodes
        self.V1, self.v2, self.NNod, self.TabRS, self.TabNTyp = self.RC.Geometry_GetNodes(
            self.RiverID, self.ReachID, self.NNod, self.TabRS, self.TabNTyp)
        self.RC.QuitRAS()

        # Set output variable IDs
        self.WSE_id = 2              # Water Surface Elevation
        self.ice_thickness_id = 184   # Ice thickness ID
        self.MinChEle = 5             # Minimum channel elevation
        self.station_nbr_id = 161     # Station number ID

        self.xs_data_modified = self.stochICE.xs_data

    def new_launch_sims(self):
        """
        Launches stochastic simulations. For each simulation, it modifies parameters, writes new files, 
        and runs HEC-RAS to compute results for water surface elevation, ice thickness, and other variables.
        """
        print('-----------------------------------------------------------')
        print('-----------------------Simulating--------------------------')
        print('-----------------------------------------------------------\n')

        times = []
        self.result_profiles = {}
        self.input_parms = pd.DataFrame()

        self.sim_keys = []
        self.flow_rates = []
        self.init_ice_thicknesses = []
        self.phi_values = []
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
            self.randomize_variables()
            
            
            self.write_new_geometry()
            
            # Store input variables
            self.sim_keys.append(self.sim_key)
            self.flow_rates.append(self.flow)
            self.init_ice_thicknesses.append(self.ice_thickness)
            self.phi_values.append(self.phi)
            self.jam_loc_upstream.append(self.location[0])
            self.jam_loc_downstream.append(self.location[1])

            # Save copies of geo and flow files
            self.store_geofile()
            self.store_flowfile()

            # Open the HEC-RAS project and run the simulation
            self.RC.Project_Open(self.stochICE.ras_file)
            self.NMsg, self.TabMsg, self.block = None, None, True

            print(f'Sim {j+1} of {self.stochICE.NSims}.')
            print(f"Running Q = {self.flow:.2f}, Phi = {self.phi:.2f}, Ice thickness = {self.ice_thickness:.2f}.")
            print(f"Ice Jam location between chainage {self.location[0]} and {self.location[1]}.")

            # Run the current simulation plan
            self.v1, self.NMsg, self.TabMsg, self.v2 = self.RC.Compute_CurrentPlan(self.NMsg, self.TabMsg, self.block)

            # Collect simulation data
            wse_list, ice_thick_list, min_channel_ele_list, stations = [], [], [], []

            # Get nodes' output data
            self.V1, self.v2, self.NNod, self.TabRS, self.TabNTyp = self.RC.Geometry_GetNodes(
                self.RiverID, self.ReachID, self.NNod, self.TabRS, self.TabNTyp)

            for i in range(self.NNod):
                if self.TabNTyp[i] == "":
                    # Water Surface Elevation
                    wse, *_ = self.RC.Output_NodeOutput(self.RiverID, self.ReachID, i+1, 0, 1, self.WSE_id)
                    wse_list.append(wse)

                    # Ice Thickness
                    thick, *_ = self.RC.Output_NodeOutput(self.RiverID, self.ReachID, i+1, 0, 1, self.ice_thickness_id)
                    ice_thick_list.append(thick)

                    # Minimum channel elevation
                    minChEle, *_ = self.RC.Output_NodeOutput(self.RiverID, self.ReachID, i+1, 0, 1, self.MinChEle)
                    min_channel_ele_list.append(minChEle)

                    # Station number
                    stationNbr, *_ = self.RC.Output_NodeOutput(self.RiverID, self.ReachID, i+1, 0, 1, self.station_nbr_id)
                    stations.append(stationNbr)

            # Store results for this simulation
            self.result_profiles[self.sim_key] = {
                'WSE': np.asarray(wse_list),
                'IceThick': np.asarray(ice_thick_list),
                'MinChEle': np.asarray(min_channel_ele_list),
                'TopIceMaxDepth': np.asarray(wse_list) - np.asarray(min_channel_ele_list),
                'station': np.asarray(stations),
            }

            # Save the 2D flood map
            tif_filename = os.path.join(
                self.stochICE.prjDir, "MonteCarlo", "SimulationTifs",
                f"WSE_{self.flow}_{self.ice_thickness}_{self.phi}_{self.porosity}_{self.location[0]}_{self.location[1]}.tif"
            )
            shutil.copyfile(self.stochICE.wse_map_path, tif_filename)

            # Remove temporary files
            for _file in glob.glob(os.path.join(self.stochICE.prjDir, "MonteCarlo", "SimulationTifs", "*.vrt")):
                os.remove(_file)

            # Stop the stopwatch and calculate simulation time
            sim_time = stopwatch.stop()
            times.append(sim_time)
            avg_time = statistics.mean(times)
            print(f'Sim time {sim_time:.2f} s, average time {avg_time:.2f} s.')
            remaining_time = (self.stochICE.NSims - (j+1)) * avg_time / 60
            print(f"Approximately {remaining_time:.2f} mins remaining in batch.\n")

        self.RC.QuitRAS()

        # Store input parameters across simulations
        self.input_parms = pd.DataFrame({
            'sim_key': self.sim_keys,
            'Q': self.flow_rates,
            'phi': self.phi_values,
            'init_ice_thicknesses': self.init_ice_thicknesses,
            'jam_loc_upstream': self.jam_loc_upstream,
            'jam_loc_downstream': self.jam_loc_downstream,
        })


    def randomize_variables(self):
        """
        Randomizes key simulation parameters such as ice thickness, phi (friction angle), 
        and ice jam location based on specified ranges in `stochICE`.
        """

        self.xs_data_modified = self.stochICE.xs_data
        
        #Flow rate
        self.flow=random.sample(self.stochICE.Q, 1)[0]
        self.set_flowrate()
       
        #phi (internal friction angle)
        self.phi=random.sample(self.stochICE.friction_angle, 1)[0]
        self.set_phi()
        
        #porosity 
        self.porosity=random.sample(self.stochICE.porosity, 1)[0]
        self.set_porosity()
        
        # thicknesses
        self.ice_thickness=random.sample(self.stochICE.Frontthick, 1)[0]
        self.thicknesses = [[f"{self.ice_thickness},{self.ice_thickness},{self.ice_thickness}"]]
        self.set_init_ice_cover_thickness()
        
        # Randomly select ice jam location from the available locations
        self.location = random.choice(self.stochICE.jam_locations)
        self.set_ice_jam_location()

    def set_init_ice_cover_thickness(self):
        """Sets the initial ice thickness for each cross-section in the geometry."""
        for item in self.xs_data_modified:
            self.xs_data_modified[item]["Ice Thickness"]['val'] = self.thicknesses[0]

    def set_phi(self):
        """Sets the ice friction angle for each cross-section in the geometry."""
        for item in self.xs_data_modified:
            self.xs_data_modified[item]['Ice Friction Angle']['val'] = self.phi

    def set_porosity(self):
        """Sets the porosity of the ice jam."""
        for item in self.xs_data_modified:
            self.xs_data_modified[item]['Ice Porosity']['val'] = self.porosity

    def set_ice_jam_location(self):
        """
        Updates the ice jam location for cross-sections where the reach matches 
        the ice jam reach and the chainage is within the specified range.
        """
        for item in self.xs_data_modified:
            if (self.xs_data_modified[item]['Reach'] == self.stochICE.ice_jam_reach and
                    self.xs_data_modified[item]['chainage'] >= self.location[0] and
                    self.xs_data_modified[item]['chainage'] <= self.location[1]):
                # Set ice jam location
                self.xs_data_modified[item]["Ice Is Channel"]['val'] = str(-1)
                self.xs_data_modified[item]["Ice Is OB"]['val'] = str(0)
            else:
                # Set default values for non-ice jam locations
                self.xs_data_modified[item]["Ice Is Channel"]['val'] = str(0)
                self.xs_data_modified[item]["Ice Is OB"]['val'] = str(0)

    def write_new_geometry(self):
        """
        Writes the new geometry with modified ice cover and ice jam locations 
        into the geometry file.
        """
        for key, item in self.xs_data_modified.items():
            try:
                self.toWrite = f"Ice Is Channel={item['Ice Is Channel']['val']}\n"
                self.lineNmb = item['Ice Is Channel']['lnNum']
                self.replace_line()
            except TypeError:
                pass

            try:
                self.toWrite = f"Ice Is OB={item['Ice Is OB']['val']}\n"
                self.lineNmb = item['Ice Is OB']['lnNum']
                self.replace_line()
            except TypeError:
                pass

            try:
                self.toWrite = f"Ice Thickness={item['Ice Thickness']['val'][0]}\n"
                self.lineNmb = item['Ice Thickness']['lnNum']
                self.replace_line()
            except TypeError:
                pass

            try:
                self.toWrite = f"Ice Porosity={item['Ice Porosity']['val'][0]}\n"
                self.lineNmb = item['Ice Porosity']['lnNum']
                self.replace_line()
            except TypeError:
                pass

            try:
                self.toWrite = f"Ice Friction Angle={item['Ice Friction Angle']['val']}\n"
                self.lineNmb = item['Ice Friction Angle']['lnNum']
                self.replace_line()
            except TypeError:
                pass

        self.close_geofile()

    def get_init_geofile_content(self):
        """Reads the initial geometry file and stores its content."""
        self.read_geofile = open(self.stochICE.geo_file, 'r')
        self.init_geofile_contents = self.read_geofile.readlines()
        self.read_geofile.close()

    def replace_line(self):
        """Replaces a line in the geometry file with modified content."""
        self.new_geofile_contents = self.init_geofile_contents
        self.new_geofile_contents[self.lineNmb] = self.toWrite

        with open(self.stochICE.geo_file, 'w') as self.out:
            self.out.writelines(self.new_geofile_contents)

    def close_geofile(self):
        """Closes the geometry file."""
        self.out.close()

    def set_flowrate(self):
        """
        Updates the flow rate in the HEC-RAS flow file. 
        This implementation assumes a single reach or a divided upstream reach.
        """
        with open(self.stochICE.flow_file, 'r') as file:
            lines = file.readlines()

        updated_lines = []
        river_reaches = []
        current_reach_index = -1

        for i, line in enumerate(lines):
            if "River Rch & RM" in line:
                current_reach_index += 1
                river_reaches.append(i)
            updated_lines.append(line)

        for idx, reach_line in enumerate(river_reaches):
            if idx == 0:
                updated_lines[reach_line + 1] = f"     {self.flow}\n"
            else:
                half_value = self.flow / 2
                updated_lines[reach_line + 1] = f"     {half_value}\n"

        with open(self.stochICE.flow_file, 'w') as file:
            file.writelines(updated_lines)

    def store_geofile(self):
        """Saves a copy of the modified geometry file with unique simulation parameters."""
        geo_copy_path = os.path.join(
            self.stochICE.prjDir, "MonteCarlo", "GeoFiles",
            f"{self.flow}_{self.ice_thickness}_{self.phi}_{self.location[0]}_{self.location[1]}.g01")
        shutil.copyfile(self.stochICE.geo_file, geo_copy_path)

    def store_flowfile(self):
        """Saves a copy of the modified flow file with unique simulation parameters."""
        flow_copy_path = os.path.join(
            self.stochICE.prjDir, "MonteCarlo", "FlowFiles",
            f"{self.flow}_{self.ice_thickness}_{self.phi}_{self.location[0]}_{self.location[1]}.f01")
        shutil.copyfile(self.stochICE.flow_file, flow_copy_path)

    def clear_results(self):
        """
        Deletes previous simulation results from the MonteCarlo folder, including *.tifs, 
        *.g0, and *.f0 files.
        """
        print('Deleting *.tifs, *.g0 and *.f0 records in MonteCarlo folder!\nFlag "clrRes=False" to suppress.\n')
        files = glob.glob(self.stochICE.tif_path + '\\*')
        for f in files:
            os.remove(f)

        files = glob.glob(self.stochICE.geo_files_path + '\\*')
        for f in files:
            os.remove(f)

        files = glob.glob(self.stochICE.flow_files_path + '\\*')
        for f in files:
            os.remove(f)

    def compress_results(self):
        """Compresses the simulation results into a zip file."""
        shutil.make_archive(
            self.stochICE.prjDir + '\\MC_results_batch_%s' % self.stochICE.ID,
            'zip',
            self.stochICE.MC_path)

    def make_ensemble_flood_map(self):
        """
        Generates an ensemble flood map by aggregating the flood maps generated from 
        multiple simulations and saves it as a TIFF file.
        """
        listdir = os.listdir(self.stochICE.tif_path)

        for count, m in enumerate(listdir):
            data_name = os.path.join(self.stochICE.tif_path, m)
            tiff = rasterio.open(data_name)
            arr = tiff.read()

            arr[arr > 0] = 1
            arr[arr < 0] = 0

            if count == 0:
                stoch = arr
            else:
                stoch = stoch + arr

        self.stoch = (stoch / self.stochICE.NSims) * 100
        self.floodmap_path = os.path.join(self.stochICE.prjDir, "FloodMaps")

        if not os.path.exists(self.floodmap_path):
            os.makedirs(self.floodmap_path)

        self.stoch = self.stoch.astype(int)

        try:
            os.remove(os.path.join(self.floodmap_path, "ensemble_floodmap.tif"))
        except FileNotFoundError:
            pass

        floodmap_output = os.path.join(self.floodmap_path, "ensemble_floodmap.tif")

        result = rasterio.open(
            floodmap_output, 'w', driver='GTiff',
            height=tiff.shape[0], width=tiff.shape[1],
            count=1, dtype=stoch.dtype, crs=tiff.crs, transform=tiff.transform)

        result.write(stoch[0, :, :], 1)
        result.close()

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
