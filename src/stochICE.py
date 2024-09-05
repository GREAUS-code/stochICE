

import os
import shutil
import numpy as np
import pickle
import time
import stochRIVICE
import stochHECRAS

class StochICE:
    """
    Class for managing stochastic simulations in ice-flood modeling 
    using both HECRAS and RIVICE models.

    Attributes:
        prjDir (str): Project directory path.
        ID (str): Batch ID for the simulations.
        ras_file (str): Path to the HEC-RAS file.
        geo_file (str): Path to the geometry file.
        flow_file (str): Path to the flow file.
        wse_map_path (str): Path to the water surface elevation map.
        NSims (int): Number of simulations.
        GrpSize (int): Group size for simulations.
        thick (list): Range of ice thickness.
        phi (list): Range of ice friction angle.
        flows (list): Range of flow rates.
        clr (bool): Whether to clear results before simulation.
        compress (bool): Whether to compress results after simulation.
        sleep (int): Sleep duration between simulations.
    """
    
    def __init__(self, prjDir: str, batch_ID: str, ras_file: str, geo_file: str, 
                 flow_file: str, wse_file: str, NSims: int, GrpSize: int, 
                 thick_range: list, phi_range: list, flow_range: list, ice_jam_reach: int,
                 ds_slope: float, max_Q: float, locations: list, 
                 riv_dwn_bc_opt: int, code: str, interval: int, days: int, 
                 timestep: int, ice_start: int, ice_end: int, profile_int: int, 
                 clrRes: bool, compRes: bool, sleep: int, stochvars: list):
        """
        Initializes the StochICE class with required parameters for stochastic simulations.

        Args:
            prjDir (str): Path to the project directory.
            batch_ID (str): Batch ID for this simulation run.
            ras_file (str): Filename of the HEC-RAS file.
            geo_file (str): Filename of the geometry file.
            flow_file (str): Filename of the flow file.
            wse_file (str): Filename for the water surface elevation map.
            NSims (int): Number of simulations to run.
            GrpSize (int): Group size for simulations.
            thick_range (list): Range of ice thickness values [min, max].
            phi_range (list): Range of ice friction angle values [min, max].
            flow_range (list): Range of flow rates [min, max].
            ds_slope (float): Downstream slope value.
            max_Q (float): Maximum discharge.
            ice_jam_reach (int): Reach with the ice jam.
            locations (list): List of ice jam locations.
            riv_dwn_bc_opt (int): Downstream boundary condition option for RIVICE.
            code (str): Simulation code, either "HECRAS" or "RIVICE".
            interval (int): Time interval for the RIVICE simulation.
            days (int): Number of days for simulation.
            timestep (int): Time step duration.
            ice_start (int): Start time for ice simulation.
            ice_end (int): End time for ice simulation.
            profile_int (int): Profile interval for simulation.
            clrRes (bool): Whether to clear previous results.
            compRes (bool): Whether to compress results after the simulation.
            sleep (int): Sleep time between simulations.
            stochvars (list): List of stochastic variables.

        Returns:
            None
        """
        self.prjDir = prjDir
        self.ID = batch_ID
        self.ras_file = os.path.join(self.prjDir, ras_file)
        self.geo_file = os.path.join(self.prjDir, geo_file)
        self.flow_file = os.path.join(self.prjDir, flow_file)
        self.geo_file_temp = os.path.join(self.prjDir, geo_file + '.temp')
        self.wse_map_path = os.path.join(self.prjDir, "MonteCarlo", wse_file)

        # Simulation parameters
        self.NSims = NSims
        self.GrpSize = GrpSize
        self.thick = thick_range
        self.phi = phi_range
        self.flows = flow_range
        self.ds_slope = ds_slope
        self.max_Q = max_Q
        self.ice_jam_reach=ice_jam_reach
        self.locations = locations
        self.riv_dwn_bc_opt = riv_dwn_bc_opt
        self.code = code
        self.interval = interval
        self.days = days
        self.timestep = timestep
        self.ice_start = ice_start
        self.ice_end = ice_end
        self.profile_int = profile_int
        self.clr = clrRes
        self.compress = compRes
        self.sleep = sleep
        self.stochvars = stochvars

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

        # HECRAS or RIVICE specific actions
        if code == 'HECRAS':
            self.stochHECRAS = stochHECRAS.StochHECRAS(self)
            self.stochHECRAS.preprocess_sims()
            self.stochHECRAS.launch_sims()

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
        print("----------------------------------------------------------")
        print(f"Running HEC-RAS project file: {self.ras_file}")
        print(f"with HEC-RAS geometry file: {self.geo_file}")
        print(f"and HEC-RAS flow file: {self.flow_file}")
        print(f"Batch ID: {self.ID}")
        print(f"Total simulations: {self.NSims}, backend: {self.code}")

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



    # def get_XS_ice_params(self):
    #     """
    #     Extracts ice parameters from cross-sections in the geometry file.

    #     Args:
    #         None

    #     Returns:
    #         None
    #     """
    #     self.xs_data = {}
    #     with open(self.geo_file, 'r') as f:
    #         flag = False
    #         for count, line in enumerate(f):
    #             if "Type RM Length L Ch R = 1" in line:
    #                 xs = line.split(",")[1].strip()
    #                 self.xs_data[xs] = {"chainage": float(xs)}
    #                 for variable in self.ice_variables:
    #                     self.xs_data[xs][variable] = {}
    #                 flag = True
    #             elif flag:
    #                 for variable in self.ice_variables:
    #                     if variable in line:
    #                         self.xs_data[xs][variable]['val'] = line.split("=")[1].strip()
    #                         self.xs_data[xs][variable]['lnNum'] = count
    #                         if len(self.xs_data[xs]) == len(self.ice_variables):
    #                             flag = False


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
                # Detect cross-section start line
                if "Type RM Length L Ch R = 1" in line:
                    xs = line.split(",")[1].strip()  # Extract cross-section identifier
                    flag = True  # Start looking for #Mann tag
    
                # Check for #Mann and prepare to read next line
                elif flag and "#Mann" in line:
                    # Read the next line to extract Manning's coefficients
                    next_line = next(f).split()  # Get the line after '#Mann'
                    print(next_line)  # Debug: Print extracted line data
                    if len(next_line) >= 8:  # Ensure there are enough numbers in the line
                        self.xs_data[xs]['Manning'] = {
                            'val_LOB': next_line[1],  # Second number
                            'val_MAIN': next_line[4],  # Fifth number
                            'val_ROB': next_line[7],  # Eighth number
                            'lnNum': count + 2  # Line number of the extracted data
                        }
                    flag = False  # Reset flag to stop looking for #Mann


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
                    xs = line.split(",")[1].strip()  # Extract cross-section identifier
                    flag = True  # Start looking for #Mann tag
    
                # Check for #Mann and prepare to read next line
                elif flag and "#Mann" in line:
                    # Read the next line to extract Bank Station coordinates
                    next_line = next(f).split()  # Get the line after '#Mann'
                    print(next_line)  # Debug: Print extracted line data
                    if len(next_line) >= 7:  # Ensure there are enough numbers in the line
                        self.xs_data[xs]['BankStations'] = {
                            'val_L': float(next_line[3]),  # Fourth number
                            'val_R': float(next_line[6]),  # Seventh number
                            'lnNum': count + 2  # Line number of the extracted data
                        }
                    flag = False  # Reset flag to stop looking for #Mann



    """
    These codes are kept because it used to work. I think that HECRAS produces a slightly different geometry file format between versions.
    """

    # def get_XS_manning(self):
    #     """
    #     Extracts Manning's roughness coefficient from the cross-sections.

    #     Args:
    #         None

    #     Returns:
    #         None
    #     """
    #     with open(self.geo_file, 'r') as f:
    #         flag = False
    #         for count, line in enumerate(f):
    #             if "Type RM Length L Ch R = 1" in line:
    #                 xs = line.split(",")[1].strip()
    #                 flag = True
    #             elif flag and "#Mann" in line:
    #                 print(line.split())
    #                 self.xs_data[xs]['Manning'] = {
    #                     'val_LOB': line.split()[1],
    #                     'val_MAIN': line.split()[4],
    #                     'val_ROB': line.split()[7],
    #                     'lnNum': count + 1
    #                 }
    #                 flag = False

    # def get_XS_bank_stations(self):
    #     """
    #     Extracts bank station coordinates from the cross-sections.

    #     Args:
    #         None

    #     Returns:
    #         None
    #     """
    #     with open(self.geo_file, 'r') as f:
    #         flag = False
    #         for count, line in enumerate(f):
    #             if "Type RM Length L Ch R = 1" in line:
    #                 xs = line.split(",")[1].strip()
    #                 flag = True
    #             elif flag and "#Mann" in line:
    #                 self.xs_data[xs]['BankStations'] = {
    #                     'val_L': float(line.split()[3]),
    #                     'val_R': float(line.split()[6]),
    #                     'lnNum': count + 1
    #                 }
    #                 flag = False

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

def save_batch(stochICE_instance: StochICE):
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

def open_batch(path: str) -> StochICE:
    """
    Loads a saved StochICE instance from a file.

    Args:
        path (str): The file path to the saved StochICE instance.

    Returns:
        StochICE: The loaded StochICE instance.
    """
    with open(path, "rb") as f:
        return pickle.load(f)

