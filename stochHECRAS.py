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


class StochHECRAS():
    
    
    def __init__(self,stochICE):
        
        self.stochICE=stochICE

        if self.stochICE.clr:
            self.clear_results()
        
    def preprocess_sims(self):

        self.RC = win32com.client.Dispatch("RAS641.HECRASCONTROLLER")
        self.RC.Project_Open(self.stochICE.ras_file)
        self.NNod, self.TabRS, self.TabNTyp = None, None, None
        self.RiverID = 1
        self.ReachID = 1
        self.V1,self.v2,self.NNod,self.TabRS,self.TabNTyp = self.RC.Geometry_GetNodes(self.RiverID, 
                                                                                      self.ReachID,
                                                                                      self.NNod,
                                                                                      self.TabRS,
                                                                                      self.TabNTyp)
        self.RC.QuitRAS()

        # HECRAS controller variable output codes (Available in Annexe E of Breaking the HEC-RAS code)
        self.WSE_id = 2
        self.ice_thickness_id = 184
        self.MinChEle= 5        


    def launch_sims(self):

        print('-----------------------------------------------------------')
        print('-----------------------Simulating--------------------------')
        print('-----------------------------------------------------------\n')
        times=[]

        self.result_profiles={}
        self.input_parms=pd.DataFrame()

        self.sim_keys=[]
        self.flow_rates=[]
        self.init_ice_thicknesses=[]
        self.phi_values=[]
        self.jam_loc_upstream=[]
        self.jam_loc_downstream=[]

        self.get_init_geofile_content()

        for j in range(0,self.stochICE.NSims):

            self.sim_key='sim_'+str(j+1)

            self.result_profiles[self.sim_key]={}

            stopwatch = Stopwatch()
            stopwatch.start()

			#modify simulation parameters
            self.randomize_variables()
            self.write_new_geometry()
            self.get_random_flowrate()

			#populate lists of input variables
            self.sim_keys.append(self.sim_key)
            self.flow_rates.append(self.flow)
            self.init_ice_thicknesses.append(self.ice_thickness)
            self.phi_values.append(self.phi)
            self.jam_loc_upstream.append(self.location[0])
            self.jam_loc_downstream.append(self.location[1])

 			#save copies of geo and flow files
            self.store_geofile()
            self.store_flowfile()

 			#HEC-RAS controller related
            self.RC.Project_Open(self.stochICE.ras_file)
            self.RC.Compute_HideComputationWindow()
            self.NMsg,self.TabMsg,self.block = None, None, True

            print('Sim %d of %d.' %(j+1, self.stochICE.NSims))
            print("Running Q = %3.2f, Phi = %3.2f, Ice thickness = %3.2f." %(self.flow,
                                                                             self.phi,
                                                                             self.ice_thickness))
            
            print("Ice Jam location between chainage %d and %d." %(self.location[0],
                                                                   self.location[1],))

            stopwatch = Stopwatch()
            stopwatch.start()
            self.v1,self.NMsg,self.TabMsg,self.v2 = self.RC.Compute_CurrentPlan(self.NMsg,
                                                                                self.TabMsg,
                                                                                self.block)
 			
 			#-----------------------
 			#get simulation data
 			#-----------------------
 			
            wse_list=[]
            ice_thick_list=[]
            min_channel_ele_list=[]

            self.V1,self.v2,self.NNod,self.TabRS,self.TabNTyp = self.RC.Geometry_GetNodes(self.RiverID,
                                                                                          self.ReachID,
                                                                                          self.NNod,
                                                                                          self.TabRS,
                                                                                          self.TabNTyp)

            # get water surface profile
            for i in range(0,self.NNod):
                if self.TabNTyp[i] =="":
                    wse,self.v1,self.v2,self.v3,self.v4,self.v5,self.v6 = self.RC.Output_NodeOutput(self.RiverID,
                                                                                                    self.ReachID,
                                                                                                    i+1,0,1,
                                                                                                    self.WSE_id)
                    wse_list.append(wse)

            # get ice thickness profile
            for i in range(0,self.NNod):
                if self.TabNTyp[i] =="":

                    thick,self.v1,self.v2,self.v3,self.v4,self.v5,self.v6 = self.RC.Output_NodeOutput(self.RiverID,
                                                                                                      self.ReachID,
                                                                                                      i+1,0,1,
                                                                                                      self.ice_thickness_id)
                    ice_thick_list.append(thick)

            # get minimum channel elevation profile
            for i in range(0,self.NNod):
                if self.TabNTyp[i] =="":
                    minChEle,self.v1,self.v2,self.v3,self.v4,self.v5,self.v6 = self.RC.Output_NodeOutput(self.RiverID,
                                                                                                         self.ReachID,
                                                                                                         i+1,0,1,
                                                                                                         self.MinChEle)
                    min_channel_ele_list.append(minChEle)

 			#-----------------------
 			#store simulation data
 			#-----------------------

            self.result_profiles[self.sim_key]['WSE']=np.asarray(wse_list)
            self.result_profiles[self.sim_key]['IceThick']=np.asarray(ice_thick_list)
            self.result_profiles[self.sim_key]['MinChEle']=np.asarray(min_channel_ele_list)
            self.result_profiles[self.sim_key]['TopIceMaxDepth']=self.result_profiles[self.sim_key]['WSE']-self.result_profiles[self.sim_key]['MinChEle']

 			#copy and store 2D flood map for 
            # tif_filename=self.stochICE.prjDir+"\\MonteCarlo\\SimulationTifs"+"\\WSE_"+str(self.flow)+"_"+str(self.ice_thickness)+"_"+str(self.phi)+".tif"
            # shutil.copyfile(self.stochICE.wse_map_path,tif_filename)
            # os.remove(self.stochICE.wse_map_path)

            vrts_to_remove = glob.glob(self.stochICE.prjDir+"\\MonteCarlo\\SimulationTifs"+"\\*.vrt")

            for _file in vrts_to_remove:

                os.remove(_file)

 			#-----------------------
 			#print sim time
 			#-----------------------

            sim_time = stopwatch.stop()
            times.append(sim_time)
            avg_time=statistics.mean(times)
            print('Sim time %3.2f s, average time %3.2f s.' %(sim_time,avg_time))
            remaining_time = (self.stochICE.NSims-(j+1))*avg_time/60
            print(f"Approximately {remaining_time:.2f} mins remaining in batch.\n")

        self.RC.QuitRAS()

		#populate input_parms with ensemble of sim input parms
        self.input_parms['sim_key']=self.sim_keys
        self.input_parms['Q']=self.flow_rates
        self.input_parms['phi']=self.phi_values
        self.input_parms['init_ice_thicknesses']=self.init_ice_thicknesses
        self.input_parms['jam_loc_upstream']=self.jam_loc_upstream
        self.input_parms['jam_loc_downstream']=self.jam_loc_downstream

        print('\nHECRAS simulation(s) complete!')
        
        
    def randomize_variables(self):

        secure_random = random.SystemRandom()

        #randomly select between min and max
        self.ice_thickness = round(secure_random.uniform(self.stochICE.thick[0], self.stochICE.thick[1]),2)
        self.phi = round(secure_random.uniform(self.stochICE.phi[0], self.stochICE.phi[1]),0)
        self.thicknesses= [[str(self.ice_thickness)+","+str(self.ice_thickness)+","+str(self.ice_thickness)]]

        self.xs_data_modified=self.stochICE.xs_data
        self.set_init_ice_cover_thickness()
        self.set_phi()

        #randomly select ice jam location from list of possible locations
        self.location = random.choice(self.stochICE.locations)
        self.set_ice_jam_location()


    def set_init_ice_cover_thickness(self):

        for item in self.xs_data_modified:
            self.xs_data_modified[item]["Ice Thickness"]['val']=self.thicknesses[0]


    def set_phi(self):

        for item in self.xs_data_modified:
            self.xs_data_modified[item]['Ice Friction Angle']['val']=self.phi


    def set_ice_jam_location(self):

        for item in self.xs_data_modified:
            if self.xs_data_modified[item]['chainage'] <= self.location[0] and self.xs_data_modified[item]['chainage'] >= self.location[1]:
                self.xs_data_modified[item]["Ice Is Channel"]['val']=str(-1)
                self.xs_data_modified[item]["Ice Is OB"]['val']=str(-1)

            else:
                self.xs_data_modified[item]["Ice Is Channel"]['val']=str(0)
                self.xs_data_modified[item]["Ice Is OB"]['val']=str(0)

        
    def write_new_geometry(self):

        for key,item in self.xs_data_modified.items():

            try:
                self.toWrite="Ice Is Channel=%s" % item['Ice Is Channel']['val']+'\n'
               
                self.lineNmb=item['Ice Is Channel']['lnNum']
                self.replace_line()

            except TypeError:
                pass

            try:
                self.toWrite="Ice Is OB=%s" % item['Ice Is OB']['val']+'\n'
                self.lineNmb=item['Ice Is OB']['lnNum']
                self.replace_line()

            except TypeError:
                pass

            try:
                self.toWrite="Ice Thickness=%s" % item['Ice Thickness']['val'][0]+'\n'
                self.lineNmb=item['Ice Thickness']['lnNum']
                self.replace_line()

            except TypeError:
                pass

            try:
                self.toWrite="Ice Friction Angle=%s" % item['Ice Friction Angle']['val']+'\n'
                self.lineNmb=item['Ice Friction Angle']['lnNum']
                self.replace_line()

            except TypeError:
                pass

       
        self.close_geofile()
       

    def get_init_geofile_content(self):
        
        self.read_geofile = open(self.stochICE.geo_file, 'r')
        self.init_geofile_contents=self.read_geofile.readlines()
        self.read_geofile.close()


    def replace_line(self):

        self.new_geofile_contents=self.init_geofile_contents
        self.new_geofile_contents[self.lineNmb] = self.toWrite

        with open(self.stochICE.geo_file, 'w') as self.out:

            self.out.writelines(self.new_geofile_contents)


    def close_geofile(self):

        self.out.close()


    def get_random_flowrate(self):
        
        secure_random = random.SystemRandom()
        self.flow = round(secure_random.uniform(self.stochICE.flows[0], self.stochICE.flows[1]),0)
        self.set_random_flowrate()


    def set_random_flowrate(self):
        """
        Places new flow rate on the sixth line of the HEC-RAS steady state flow file
        """

        with open(self.stochICE.flow_file,'r') as txt:
            text=txt.readlines()
            text[5]='     %s\n'%self.flow

            with open(self.stochICE.flow_file,'w') as txt:
                txt.writelines(text)

                
    def store_geofile(self):

        geo_copy_path=self.stochICE.prjDir+"\\MonteCarlo\\GeoFiles\\"+str(self.flow)+"_"+str(self.ice_thickness)+"_"+str(self.phi)+"_"+str(self.location[0])+"_"+str(self.location[1])+".g01"
        shutil.copyfile(self.stochICE.geo_file,geo_copy_path)


    def store_flowfile(self):

        flow_copy_path=self.stochICE.prjDir+"\\MonteCarlo\\FlowFiles\\"+str(self.flow)+"_"+str(self.ice_thickness)+"_"+str(self.phi)+"_"+str(self.location[0])+"_"+str(self.location[1])+".f01"
        shutil.copyfile(self.stochICE.flow_file,flow_copy_path)


    def clear_results(self):

        print('Deleting *.tifs, *.g0 and *.f0 records in MonteCarlo folder!\nFlag "clrRes=False" to suppress.\n')
        files = glob.glob(self.stochICE.tif_path+'\\*')
        for f in files:
            os.remove(f)

        files = glob.glob(self.stochICE.geo_files_path+'\\*')
        for f in files:
            os.remove(f)

        files = glob.glob(self.stochICE.flow_files_path+'\\*')
        for f in files:
            os.remove(f)


    def compress_results(self):

        shutil.make_archive(self.stochICE.prjDir+'\\MC_results_batch_%s'%self.stochICE.ID, 
                            'zip', 
                            self.stochICE.MC_path)


    def make_ensemble_flood_map(self):

        listdir = os.listdir(self.stochICE.tif_path)

        for count,m in enumerate(listdir):

            data_name = self.stochICE.tif_path + "\\"+m
            tiff = rasterio.open(data_name)
            arr = tiff.read()

            arr[arr > 0] = 1
            arr[arr < 0] = 0

            #Check if first image
            if count == 0:
                stoch=arr
            else:
                stoch = stoch + arr

        self.stoch=(stoch/self.stochICE.NSims)*100
        self.floodmap_path=self.stochICE.prjDir+"\FloodMaps"

        if not os.path.exists(self.floodmap_path):
            os.makedirs(self.floodmap_path)

        self.stoch = self.stoch.astype(int)

        try:
            os.remove(self.floodmap_path+"\ensemble_floodmap.tif")#delete previous map
        except FileNotFoundError:
            pass

        floodmap_output = self.floodmap_path+"\ensemble_floodmap.tif"
        
        result = rasterio.open(floodmap_output,'w',driver='GTiff',
                               height=tiff.shape[0],
                               width=tiff.shape[1],
                               count=1,
                               dtype=stoch.dtype,
                               crs=tiff.crs,
                               transform=tiff.transform)
        
        result.write(stoch[0,:,:],1)
        result.close()






class Stopwatch:
    def __init__(self):
        self.start_time = None

    def start(self):
        self.start_time = time.time()

    def stop(self):
        if self.start_time is None:
            raise ValueError("Stopwatch has not been started.")
        elapsed_time = time.time() - self.start_time
        self.start_time = None
        return elapsed_time
