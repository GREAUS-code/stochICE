# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 09:26:34 2023

@author: Jason
"""


import os
import glob
import shutil
import random
import time
import win32com.client
import numpy as np
import statistics
#import rasterio
import pandas as pd
import re


import stochRIVICE

# 
class stochICE():

    def __init__(self,prjDir,
                      batch_ID,
                      ras_file,
                      geo_file,
                      flow_file,
                      wse_file,
                      NSims,
                      thick_range,
                      phi_range,
                      flow_range,
                      locations,
                      code,
                      clrRes,
                      compRes):

        #paths
        self.prjDir = prjDir
        self.ID=batch_ID
        self.ras_file=self.prjDir +'\\'+ ras_file
        self.geo_file=self.prjDir +'\\'+ geo_file
        self.geo_file_temp=self.prjDir +'\\'+ geo_file + '.temp'
        self.flow_file=self.prjDir +'\\'+ flow_file
        self.wse_map_path=self.prjDir + '\\MonteCarlo\\' + wse_file

        #variables
        self.NSims=NSims
        self.thick=thick_range
        self.phi=phi_range
        self.flows=flow_range
        self.locations=locations
        self.code=code
        self.clr=clrRes
        self.compress=compRes
        

        #ice variables in geo file
        self.iceVariables=['Ice Thickness',
                   'Ice Mann',
                   'Ice Specific Gravity',
                   'Ice Is Channel',
                   'Ice Is OB',
                   'Ice Friction Angle',
                   'Ice Porosity',
                   'Ice K1',
                   'Ice Max Mean',
                   'Ice Cohesion',
                   'Ice Fixed Mann',
                   ]


		#Add functions here to run them by default

        self.printHeader()
        self.setupMonteCarloDir()

        # if self.clr:
        #     self.clearResults()

        # Common functions for either HECRAS or RIVICE
        self.preprocessSimulations()
        self.getXSectionIceParameters()
        self.getXSectionManning()
        self.getXSectionBankStations()
        self.getXSectionGeometry()
        self.getXSectionMainChannelGeometry()
        self.getBridgeData()


        """
        HECRAS specific functions will go here
        """

        if code =='HECRAS':
            print('Using HECRAS backend')
            # self.launch_HECRAS_simulations()

        """
        RIVICE specific functions will go here
        """

        if code =='RIVICE':
        
            print('Using RIVICE backend')
            
            #create stochRIVICE instance
            self.stochRIVICE=stochRIVICE.StochRIVICE(self,5)

            self.stochRIVICE.get_riv_xsections()
            self.stochRIVICE.RiviceChainage()
            self.stochRIVICE.RiviceManning()
            self.stochRIVICE.RiviceXSectionMainChannelGeometry()
            self.stochRIVICE.RiviceReach()
            self.stochRIVICE.RiviceXSectionID()
            self.stochRIVICE.RiviceDistanceXSectionPrecedente()
            self.stochRIVICE.WriteRiviceCd1test()



        # if self.compress:
        #     self.compressResults()


    def printHeader(self):


        print("\n")
        print("----------------------------------------------------------")
        print("Running HEC-RAS iceJammer2000 flood modelling system v0.1")
        print("----------------------------------------------------------")

        print(r"""⠀⠀⠀
            ⣀⣤⣴⣶⣶⣶⣦⣤⣀⠀⠀⠀⣠⣴⣿⣿⣿⣿⣿⣿⣿⣿⣶⣄⡀⠀⠀
        ⠀⣴⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠇⣠⣾⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣦⠀
        ⠀⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠏⢰⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠀
        ⠀⣿⣿⠿⠿⠿⠿⠟⠛⠛⠛⠀⠛⠛⠛⠛⠛⠛⠛⠻⠿⠿⠿⠿⠿⠿⠿⠿⠿⠀
        ⠀⠀⠀⠀⠀⠀⠀⢠⣤⣶⣶⣶⠋⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
        ⠀⠀⠀⠀⠀⠀⠀⠀⠙⠻⢿⣿⣿⣶⣶⣦⣤⣄⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
        ⠀⠀⠀⠀⠀⠀⣀⣀⣀⣤⣤⣤⣽⣿⣿⣿⣿⡿⠇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
        ⠀⠀⣤⣶⣿⣿⣿⣿⣿⣿⣿⡿⠛⠛⠋⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
        ⠀⠀⠙⠛⠿⢿⣯⣭⣝⡛⠻⢿⣿⣿⣷⣶⣶⣦⣤⣤⣄⡀⠀⠀⠀⠀⠀⠀⠀⠀
        ⠀⠀⠀⠀⠀⠀⠀⠉⠙⠛⠿⢶⣾⣿⣿⣿⣿⣿⣿⣿⡿⣿⣷⠀⠀⠀⠀⠀⠀⠀
        ⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣀⣤⣾⣿⣿⢿⣿⣿⢟⣡⣼⣿⠟⠀⠀⠀⠀⠀⠀⠀
        ⠀⠀⠀⠀⠀⠀⠀⣠⣶⣿⣿⣿⠟⣡⣾⣿⣿⣿⣿⡿⠋⠁⠀⠀⠀⠀⠀⠀⠀⠀
        ⠀⠀⠀⠀⠀⠀⢰⣿⣿⠹⣿⣿⣄⠻⣿⣿⣿⠻⣿⣿⣦⣄⣀⠀⠀⠀⠀⠀⠀⠀
        ⠀⠀⠀⠀⠀⠀⠘⠛⠛⠓⠈⠛⠛⠛⠊⠛⠛⠓⠀⠙⠛⠛⠛⠛⠓⠒""")


        print("Running HEC-RAS project file:")
        print("%s\n" % self.ras_file)
        print("with HEC-RAS geometry file:")
        print("%s\n" % self.geo_file)
        print("and HEC-RAS flow file:")
        print("%s\n" % self.flow_file)

        print('                 ------------------\n')
        print("                 Batch ID is: %s\n" % self.ID)
        print('                 ------------------\n')

        print("A total of %s simulations will be run." % str(self.NSims))
        print("Ice thickness varies between %3.2f and %3.2f m" %(self.thick[0],self.thick[1]))
        print("Discharges vary between %3.2f and %3.2f m**3/s" %(self.flows[0],self.flows[1]))
        print("Phi varies between %3.2f and %3.2f\n" %(self.phi[0],self.phi[1]))
        print("Lodgement location randomly chosen between %2.0f provided locations.\n" %(len(self.locations)))


        print("""
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@%...    .  ./&@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@@@@@@,..   @@@@@@/.  #@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@@@/   ..,@@&  @@@@@@ . /@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@      /@@@&(#@@   #@  .  @@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@, .   . @@&@@*  &%/ ..     &@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@/   .   (@@/@@(, @(@@ .     &@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@%        @@@@.(@.  . ...      @@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@*..       .. .%@,.           . @@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@.        ,.   . .   .      @    @@@@@@@@@@@@@@@
@@@@@@@@@@@@@@(@.      .(.  ..   (@    . *&    @@@@@@@@@@@@@@@
@@@@@@@@@@@@@#&             .,%, .      @ .  ( @@@@@@@@@@@@@@@
@@@@@@@@@@@@@ .     .  @    @ .       # . .  @@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@* .   .  %    @@/        # @@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@,      .@@& .@&. ..      & @@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@        @@@@%%         .  %, @@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@*       ,@@@@& .&%. (       &@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@       &@(@@@@@/. @&.  %@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@      @@.@@@@@*#. (&..*@@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@.   &%@@&#@@@@/%(  /@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@(. %@@@ @@@@@@@@  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@*@@@@@@@@@@@@@% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

                    You ready to jam?
			                                                                    """)

    def setupMonteCarloDir(self):

        self.MC_path=self.prjDir+"\MonteCarlo"
        self.tif_path=self.MC_path+"\SimulationTifs"
        self.geoFiles_path=self.MC_path+"\GeoFiles"
        self.flowFiles_path=self.MC_path+"\FlowFiles"

        shutil.copyfile(self.geo_file,self.geo_file_temp)

        if not os.path.exists(self.MC_path):
            os.makedirs(self.MC_path)
            print("Created .\MonteCarlo")

        if not os.path.exists(self.tif_path):
            os.makedirs(self.tif_path)
            print("Created .\MonteCarlo\SimulationTifs")

        if not os.path.exists(self.geoFiles_path):
            os.makedirs(self.geoFiles_path)
            print("Created .\MonteCarlo\GeoFiles")

        if not os.path.exists(self.flowFiles_path):
            os.makedirs(self.flowFiles_path)
            print("Created .\MonteCarlo\FlowFiles")


    def preprocessSimulations(self):

        #Import HEC-RAS Controller
        self.RC = win32com.client.Dispatch("RAS641.HECRASCONTROLLER")
        # self.RC = win32com.client.Dispatch("RAS631.HECRASCONTROLLER")

        self.RC.Project_Open(self.ras_file)
        self.NNod, self.TabRS, self.TabNTyp = None, None, None
        self.RiverID = 1
        self.ReachID = 1
        self.V1,self.v2,self.NNod,self.TabRS,self.TabNTyp = self.RC.Geometry_GetNodes(self.RiverID, self.ReachID,self.NNod,self.TabRS,self.TabNTyp)
        self.RC.QuitRAS()

        # HECRAS controller variable output codes (Available in Annexe E of Breaking the HEC-RAS code)
        self.WSE_id = 2
        self.ice_thick_id = 184
        self.MinChEle= 5



    def clearResults(self):

        print('Deleting *.tifs, *.g0 and *.f0 records in MonteCarlo folder!\nFlag "clrRes=False" to suppress.\n')
        files = glob.glob(self.tif_path+'\\*')
        for f in files:
            os.remove(f)

        files = glob.glob(self.geoFiles_path+'\\*')
        for f in files:
            os.remove(f)

        files = glob.glob(self.flowFiles_path+'\\*')
        for f in files:
            os.remove(f)


    def compressResults(self):

        shutil.make_archive(self.prjDir+'\\MC_results_batch_%s'%self.ID, 'zip', self.MC_path)


    def launch_HECRAS_simulations(self):

        print('-----------------------------------------------------------')
        print('-----------------------Simulating--------------------------')
        print('-----------------------------------------------------------\n')
        times=[]

        self.simProfileData={}
        self.simInputPars=pd.DataFrame()

        self.simKeys=[]
        self.flowRates=[]
        self.init_thick=[]
        self.phi_values=[]
        self.locUp=[]
        self.locDw=[]

        self.getInitGeoFileContents()

        for j in range(0,self.NSims):

            self.simKey='sim_'+str(j+1)

            self.simProfileData[self.simKey]={}

            stopwatch = Stopwatch()
            stopwatch.start()

			#modify simulation parameters
            self.randomlyModifyGeometryFile()
            self.writeNewGeometry()
            self.randomlyModifyFlowFile()

			#populate lists of input variables
            self.simKeys.append(self.simKey)
            self.flowRates.append(self.flow)
            self.init_thick.append(self.ice_thick)
            self.phi_values.append(self.PHI)
            self.locUp.append(self.location[0])
            self.locDw.append(self.location[1])

			#save copies of geo and flow files
            self.storeGeoFile()
            self.storeFlowFile()

			#HEC-RAS controller related
            self.RC.Project_Open(self.ras_file)
            self.RC.Compute_HideComputationWindow()
            self.NMsg,self.TabMsg,self.block = None, None, True

            print('Sim %d of %d.' %(j+1, self.NSims))
            print("Running Q = %3.2f, Phi = %3.2f, Ice thickness = %3.2f." %(self.flow,self.PHI,self.ice_thick))
            print("Ice Jam location between chainage %d and %d." %(self.location[0],self.location[1],))

            stopwatch = Stopwatch()
            stopwatch.start()
            self.v1,self.NMsg,self.TabMsg,self.v2 = self.RC.Compute_CurrentPlan(self.NMsg,self.TabMsg,self.block)
			
            
			#-----------------------
			#get simulation data
			#-----------------------
			
            wseList=[]
            iceThickList=[]
            minChEleList=[]

            self.V1,self.v2,self.NNod,self.TabRS,self.TabNTyp = self.RC.Geometry_GetNodes(self.RiverID, self.ReachID,self.NNod,self.TabRS,self.TabNTyp)

            # get water surface profile
            for i in range(0,self.NNod):
                if self.TabNTyp[i] =="":
                    wse,self.v1,self.v2,self.v3,self.v4,self.v5,self.v6 = self.RC.Output_NodeOutput(self.RiverID, self.ReachID,i+1,0,1,self.WSE_id)
                    wseList.append(wse)

            # get ice thickness profile
            for i in range(0,self.NNod):
                if self.TabNTyp[i] =="":

                    thick,self.v1,self.v2,self.v3,self.v4,self.v5,self.v6 = self.RC.Output_NodeOutput(self.RiverID, self.ReachID,i+1,0,1,self.ice_thick_id)
                    iceThickList.append(thick)

            # get minimum channel elevation profile
            for i in range(0,self.NNod):
                if self.TabNTyp[i] =="":
                    minChEle,self.v1,self.v2,self.v3,self.v4,self.v5,self.v6 = self.RC.Output_NodeOutput(self.RiverID, self.ReachID,i+1,0,1,self.MinChEle)
                    minChEleList.append(minChEle)



			#-----------------------
			#store simulation data
			#-----------------------

            self.simProfileData[self.simKey]['WSE']=np.asarray(wseList)
            self.simProfileData[self.simKey]['IceThick']=np.asarray(iceThickList)
            self.simProfileData[self.simKey]['MinChEle']=np.asarray(minChEleList)
            self.simProfileData[self.simKey]['TopIceMaxDepth']=self.simProfileData[self.simKey]['WSE']-self.simProfileData[self.simKey]['MinChEle']

			#copy and store 2D flood map for 
            # tif_filename=self.prjDir+"\\MonteCarlo\\SimulationTifs"+"\\WSE_"+str(self.flow)+"_"+str(self.ice_thick)+"_"+str(self.PHI)+".tif"
            # shutil.copyfile(self.wse_map_path,tif_filename)
            # os.remove(self.wse_map_path)

            # vrts_to_remove = glob.glob(self.prjDir+"\\MonteCarlo\\SimulationTifs"+"\\*.vrt")

            # for _file in vrts_to_remove:

            #     os.remove(_file)


			#-----------------------
			#print sim time
			#-----------------------

            sim_time = stopwatch.stop()
            times.append(sim_time)
            avgTime=statistics.mean(times)
            print('Sim time %3.2f s, average time %3.2f s.' %(sim_time,avgTime))
            remainingTime = (self.NSims-(j+1))*avgTime/60
            print(f"Approximately {remainingTime:.2f} mins remaining in batch.\n")

        self.RC.QuitRAS()

		#populate simInputPars with ensemble of sim input parms
        self.simInputPars['simKey']=self.simKeys
        self.simInputPars['Q']=self.flowRates
        self.simInputPars['phi']=self.phi_values
        self.simInputPars['init_thick']=self.init_thick
        self.simInputPars['locUp']=self.locUp
        self.simInputPars['locDw']=self.locDw

        print('\nSimulations complete!')

    def getXSectionIceParameters(self):

        self.xsData={}

        with open(self.geo_file, 'r') as f:
			
            flag=False

            for count,line in enumerate(f):

                if "Type RM Length L Ch R = 1" in line:
					
                    xs=line.split(",")[1].strip()
                    self.xsData[xs]={}
                    self.xsData[xs]['chainage']=float(xs)
                    
                    for variable in self.iceVariables:
                        self.xsData[xs][variable]={}
                    
                    flag=True
                    variableCounter=0
					
                else:
                    
                    if flag:
                        
                        for count, variable in enumerate(self.iceVariables):
                            
                            if variable in line:
                                
                                variableCounter = variableCounter + 1
                                
                                self.xsData[xs][variable]['val']=line.split("=")[1].strip()
                                self.xsData[xs][variable]['lnNum']=count
                                
                            if variableCounter == len(self.iceVariables):
                                flag=False
                                                         

    def getXSectionManning(self):

        with open(self.geo_file, 'r') as f:
			
            flag, flag2 = False, False

            for count,line in enumerate(f):

                if "Type RM Length L Ch R = 1" in line:

                    xs=line.split(",")[1].strip()
                    
                    flag=True
                    
                elif flag2:
                    
                    self.xsData[xs]['Manning'] = {}
                    self.xsData[xs]['Manning']['val_LOB'] = line.split()[1]
                    self.xsData[xs]['Manning']['val_MAIN'] = line.split()[4]
                    self.xsData[xs]['Manning']['val_ROB'] = line.split()[7]
                    self.xsData[xs]['Manning']['lnNum'] = count + 1
                    
                    flag, flag2 = False, False            
                    
                else:
                    
                    if flag:
                            
                        if '#Mann' in line:
                            
                            flag2 = True
                            
                            
                            
    def getXSectionBankStations(self):

        with open(self.geo_file, 'r') as f:
			
            flag, flag2 = False, False

            for count,line in enumerate(f):

                if "Type RM Length L Ch R = 1" in line:

                    xs=line.split(",")[1].strip()
                    
                    flag=True
                    
                elif flag2:
                    
                    self.xsData[xs]['BankStations'] = {}
                    self.xsData[xs]['BankStations']['val_L'] = float(line.split()[3])
                    self.xsData[xs]['BankStations']['val_R'] = float(line.split()[6])
                    self.xsData[xs]['BankStations']['lnNum'] = count + 1
                    
                    flag, flag2 = False, False            
                    
                else:
                    
                    if flag:
                            
                        if '#Mann' in line:
                            
                            flag2 = True
                            
                            
    def getXSectionGeometry(self):
    
        with open(self.geo_file, 'r') as f:
    			
            flag, flag2 = False, False
    
            for count,line in enumerate(f):
    
                if "Type RM Length L Ch R = 1" in line:
    
                    xs=line.split(",")[1].strip()
                    
                    self.xsData[xs]['Geometry'] = {}
                    
                    flag=True
                    
                    
                elif '#Mann' in line:
                    
                    self.xsData[xs]['Geometry']['lnNum_end'] = count
                    
                    self.xsData[xs]['Geometry']['xy'] = np.float_(self.xsData[xs]['Geometry']['xy'])
                    
                    flag, flag2 = False, False
                    
                   
                elif flag2:
                    
                    xy = []
                    
                    for i in range(10):
                        
                        if line[8*i:8*i+8] == '\n':
                        
                            break
                        
                        xy.append(line[8*i:8*i+8])
                        
                    
                    self.xsData[xs]['Geometry']['xy'] = self.xsData[xs]['Geometry']['xy'] + xy
                             
                    
                else:
                    
                    if flag:
                            
                        if '#Sta/Elev' in line:
                            
                            self.xsData[xs]['Geometry']['xy'] = []
                            self.xsData[xs]['Geometry']['lnNum_start'] = count + 2
                            
                            flag2 = True
                            
                            
    def getXSectionMainChannelGeometry(self):
        
        self.getXSectionGeometry()
        
        for xs in self.xsData:       
            
            self.xsData[xs]['MainChannelGeometry'] = {}
            self.xsData[xs]['MainChannelGeometry']['xy'] = []
            
            for i in range(0,len(self.xsData[xs]['Geometry']['xy']),2):
                
                if self.xsData[xs]['Geometry']['xy'][i] >= self.xsData[xs]['BankStations']['val_L'] and self.xsData[xs]['Geometry']['xy'][i] <= self.xsData[xs]['BankStations']['val_R'] : 
                    
                    self.xsData[xs]['MainChannelGeometry']['xy'].append(self.xsData[xs]['Geometry']['xy'][i])
                    self.xsData[xs]['MainChannelGeometry']['xy'].append(self.xsData[xs]['Geometry']['xy'][i+1])

    
    def getBridgeData(self):
        
        self.Bridge = False
        
        self.BridgeData = {}
        
        Bridge_number = 1
        
        with open(self.geo_file, 'r') as f:    			
    
            for count,line in enumerate(f):
    
                if "Type RM Length L Ch R = 3" in line:
                    
                    self.Bridge = True
    
                    xs=line.split(",")[1].strip()
                    self.BridgeData[str(Bridge_number)]={}  
                    self.BridgeData[str(Bridge_number)]['chainage'] = float(xs)
                    
                    Bridge_number += 1

                    


    def randomlyModifyGeometryFile(self):

        secure_random = random.SystemRandom()

        #randomly select between min and max
        self.ice_thick = round(secure_random.uniform(self.thick[0], self.thick[1]),2)
        self.PHI = round(secure_random.uniform(self.phi[0], self.phi[1]),0)
        self.thicknesses= [[str(self.ice_thick)+","+str(self.ice_thick)+","+str(self.ice_thick)]]

        self.xsData_mod=self.xsData
        self.modifyIceCoverThickness()
        self.modifyPhi()

        #randomly select ice jam location from list of possible locations
        self.location = random.choice(self.locations)
        self.modifyIceJamLocation()


    def modifyIceCoverThickness(self):

        for item in self.xsData_mod:
            self.xsData_mod[item]["Ice Thickness"]['val']=self.thicknesses[0]


#               for count, reachBounds in enumerate(bounds):
#
#                   if xsData[item]['chainage'] <= reachBounds[0] and xsData[item]['chainage'] >= reachBounds[1]:
#                       xsData[item]["Ice Thickness"]['val']=thicknesses[count][0]


    def modifyPhi(self):

        for item in self.xsData_mod:
            self.xsData_mod[item]['Ice Friction Angle']['val']=self.PHI

#               for count, reachBounds in enumerate(bounds):
#
#                   if xsData[item]['chainage'] <= reachBounds[0] and xsData[item]['chainage'] >= reachBounds[1]:
#                       xsData[item]['Ice Friction Angle']['val']=phi[count]


    def modifyIceJamLocation(self):

        for item in self.xsData_mod:
            if self.xsData_mod[item]['chainage'] <= self.location[0] and self.xsData_mod[item]['chainage'] >= self.location[1]:
                self.xsData_mod[item]["Ice Is Channel"]['val']=str(-1)
                self.xsData_mod[item]["Ice Is OB"]['val']=str(-1)

            else:
                self.xsData_mod[item]["Ice Is Channel"]['val']=str(0)
                self.xsData_mod[item]["Ice Is OB"]['val']=str(0)


    def writeNewGeometry(self):

        for key,item in self.xsData_mod.items():

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

       
        self.closeGeoFile()
       

    def getInitGeoFileContents(self):
        self.readGeoFile = open(self.geo_file, 'r')
        self.initGeoFileContents=self.readGeoFile.readlines()
        self.readGeoFile.close()


    def replace_line(self):

        self.newGeoFileContents=self.initGeoFileContents
        self.newGeoFileContents[self.lineNmb] = self.toWrite

        with open(self.geo_file, 'w') as self.out:

            self.out.writelines(self.newGeoFileContents)

    def closeGeoFile(self):

        self.out.close()

    def randomlyModifyFlowFile(self):
        secure_random = random.SystemRandom()
        self.flow = round(secure_random.uniform(self.flows[0], self.flows[1]),0)
        self.replaceFlow()

    def replaceFlow(self):
        """
        Places new flow rate on the sixth line of the HEC-RAS steady state flow file
        """

        with open(self.flow_file,'r') as txt:
            text=txt.readlines()
            text[5]='     %s\n'%self.flow

            with open(self.flow_file,'w') as txt:
                txt.writelines(text)

    def storeGeoFile(self):

        geoCopyPath=self.prjDir+"\\MonteCarlo\\GeoFiles\\"+str(self.flow)+"_"+str(self.ice_thick)+"_"+str(self.PHI)+"_"+str(self.location[0])+"_"+str(self.location[1])+".g01"
        shutil.copyfile(self.geo_file,geoCopyPath)

    def storeFlowFile(self):

        flowCopyPath=self.prjDir+"\\MonteCarlo\\FlowFiles\\"+str(self.flow)+"_"+str(self.ice_thick)+"_"+str(self.PHI)+"_"+str(self.location[0])+"_"+str(self.location[1])+".f01"
        shutil.copyfile(self.flow_file,flowCopyPath)

    def createEnsembleFloodMap(self):

        listdir = os.listdir(self.tif_path)

        for count,m in enumerate(listdir):

            data_name = self.tif_path + "\\"+m
            tiff = rasterio.open(data_name)
            arr = tiff.read()

            arr[arr > 0] = 1
            arr[arr < 0] = 0

            #Check if first image
            if count == 0:
                stoch=arr
            else:
                stoch = stoch + arr

        # stoch=(stoch/NSimul)*100
        self.stoch=(stoch/self.NSims)*100


        self.floodmap_path=self.prjDir+"\FloodMaps"


        if not os.path.exists(self.floodmap_path):
            os.makedirs(self.floodmap_path)

        self.stoch = self.stoch.astype(int)

        try:

            os.remove(self.floodmap_path+"\ensemble_floodmap.tif")#delete previous map

        except FileNotFoundError:
            pass

        floodmap_output = self.floodmap_path+"\ensemble_floodmap.tif"
        result = rasterio.open(floodmap_output,'w',driver='GTiff',height=tiff.shape[0],width=tiff.shape[1],count=1,dtype=stoch.dtype,crs=tiff.crs,transform=tiff.transform)
        result.write(stoch[0,:,:],1)
        result.close()
        
    #______________________________FONCTIONS DE CONVERSION HEC-RAS VERS RIVICE______________________________#
        
    
    # def RiviceParameters(self):
        
    #     self.RiviceParameters = {}
    #     self.RiviceParameters['Rivice XS Interpolation Interval'] = 5
        
        
     
    # def RiviceXSections(self): 
        
    #     self.RivicexsData = {}
        
    #     xs_number = 1
        
    #     xs_precedente = ''
        
    #     for xs in self.xsData:
            
    #         if xs_precedente == '':
            
    #             self.RivicexsData[str(xs_number)] = {}
    #             self.RivicexsData[str(xs_number)]['Hecras xs'] = xs
                       
    #         elif self.xsData[xs]['Manning']['val_MAIN'] != self.xsData[xs_precedente]['Manning']['val_MAIN'] :
                
    #             self.RivicexsData[str(xs_number)] = {}
    #             self.RivicexsData[str(xs_number)]['Hecras xs'] = xs_precedente
                
    #             xs_number += 1
                
    #             self.RivicexsData[str(xs_number)] = {}
    #             self.RivicexsData[str(xs_number)]['Hecras xs'] = xs
            
    #         else :
                
    #             self.RivicexsData[str(xs_number)] = {}
    #             self.RivicexsData[str(xs_number)]['Hecras xs'] = xs
        
    #         xs_number += 1        
    #         xs_precedente = xs
            
            
    # def RiviceChainage(self):
        
    #     def mon_arrondi(x, base):
    #         return base * round(x/base)
        
    #     flag = False
    #     reste_arrondi = 0            
    #     xs_number_precedent = ''
    #     xs_number_precedent_precedent = ''
        
    #     for xs_number in self.RivicexsData:
            
    #         if xs_number == '1':
    #             self.RivicexsData[xs_number]['Rivice Chainage'] = 0
                
    #         elif xs_number == str(len(self.RivicexsData)):
    #             self.RivicexsData[xs_number]['Rivice Chainage'] = mon_arrondi(self.RivicexsData[xs_number_precedent]['Rivice Chainage'] + (int(self.RivicexsData[xs_number_precedent]['Hecras xs']) - int(self.RivicexsData[xs_number]['Hecras xs'])), self.RiviceParameters['Rivice XS Interpolation Interval'])
                 
    #         elif flag:
    #             self.RivicexsData[xs_number]['Rivice Chainage'] = self.RivicexsData[xs_number_precedent]['Rivice Chainage'] + (int(self.RivicexsData[xs_number_precedent]['Hecras xs']) - (reste_arrondi + int(self.RivicexsData[xs_number]['Hecras xs'])))
    #             flag = False
                
    #         elif self.RivicexsData[xs_number]['Hecras xs'] == self.RivicexsData[xs_number_precedent]['Hecras xs'] :
                
    #             self.RivicexsData[xs_number]['Rivice Chainage'] = mon_arrondi(self.RivicexsData[xs_number_precedent_precedent]['Rivice Chainage'] + (int(self.RivicexsData[xs_number_precedent_precedent]['Hecras xs']) - int(self.RivicexsData[xs_number_precedent]['Hecras xs'])), self.RiviceParameters['Rivice XS Interpolation Interval'])
                
    #             self.RivicexsData[xs_number_precedent]['Rivice Chainage'] = self.RivicexsData[xs_number]['Rivice Chainage']   
                
    #             reste_arrondi = mon_arrondi(self.RivicexsData[xs_number_precedent_precedent]['Rivice Chainage'] + (int(self.RivicexsData[xs_number_precedent_precedent]['Hecras xs']) - int(self.RivicexsData[xs_number_precedent]['Hecras xs'])), self.RiviceParameters['Rivice XS Interpolation Interval']) - (self.RivicexsData[xs_number_precedent_precedent]['Rivice Chainage'] + (int(self.RivicexsData[xs_number_precedent_precedent]['Hecras xs']) - int(self.RivicexsData[xs_number_precedent]['Hecras xs'])))
    #             flag = True                                                                                    
            
    #         else:
                
    #             self.RivicexsData[xs_number]['Rivice Chainage'] = self.RivicexsData[xs_number_precedent]['Rivice Chainage'] + (int(self.RivicexsData[xs_number_precedent]['Hecras xs']) - int(self.RivicexsData[xs_number]['Hecras xs']))        
            
    #         xs_number_precedent_precedent = xs_number_precedent
            
    #         xs_number_precedent = xs_number
    
     
    # def RiviceManning(self):
        
    #     flag = False
        
    #     xs_number_precedent = ''
        
    #     for xs_number in self.RivicexsData:
            
            
    #         if xs_number == '1':
                
    #             self.RivicexsData[xs_number]['Manning'] = self.xsData[self.RivicexsData[xs_number]['Hecras xs']]['Manning']['val_MAIN']
            
            
    #         elif self.RivicexsData[xs_number]['Hecras xs'] == self.RivicexsData[xs_number_precedent]['Hecras xs']:
                
    #             flag = True
                
                
    #         elif flag:
                
    #             self.RivicexsData[xs_number_precedent]['Manning'] = self.xsData[self.RivicexsData[xs_number]['Hecras xs']]['Manning']['val_MAIN']
    #             self.RivicexsData[xs_number]['Manning'] = self.xsData[self.RivicexsData[xs_number]['Hecras xs']]['Manning']['val_MAIN']
                
    #             flag = False
                
                
    #         else:
            
    #             self.RivicexsData[xs_number]['Manning'] = self.xsData[self.RivicexsData[xs_number]['Hecras xs']]['Manning']['val_MAIN']


    #         xs_number_precedent = xs_number
            
            
            
    # def RiviceXSectionMainChannelGeometry(self):
        
    #     for xs_number in self.RivicexsData:
            
    #         self.RivicexsData[xs_number]['MainChannelGeometry'] = self.xsData[self.RivicexsData[xs_number]['Hecras xs']]['MainChannelGeometry']['xy']
            
            
        
    # def RiviceReach(self):
        
    #     Reach = 1
        
    #     xs_number_precedent = ''
        
    #     for xs_number in self.RivicexsData:
            
    #         if xs_number == '1':
            
    #             self.RivicexsData[xs_number]['Reach'] = Reach
                
                
    #         elif self.RivicexsData[xs_number]['Manning'] != self.RivicexsData[xs_number_precedent]['Manning']:
                
    #             Reach +=1
                
    #             self.RivicexsData[xs_number]['Reach'] = Reach
                
    #         else:
                
    #             self.RivicexsData[xs_number]['Reach'] = Reach
        
            
    #         xs_number_precedent = xs_number
            
            
    # def RiviceXSectionID(self):
        
    #     for xs_number in self.RivicexsData:
            
    #         self.RivicexsData[xs_number]['XS ID'] = int(xs_number)*1000
        
        
        
        
    # def RiviceDistanceXSectionPrecedente(self):
        
    #     xs_number_precedent = ''
        
    #     for xs_number in self.RivicexsData:
            
    #         if xs_number == '1':
            
    #             self.RivicexsData[xs_number]['Distance XS Precedente'] = 0
            
    #         else:
                
    #             self.RivicexsData[xs_number]['Distance XS Precedente'] = self.RivicexsData[xs_number]['Rivice Chainage'] - self.RivicexsData[xs_number_precedent]['Rivice Chainage']
                
                
    #         xs_number_precedent =xs_number   
                
                
   

    # def WriteRiviceCd1test(self):    
        
        
    #     Cd1test = open(self.prjDir + '/' + 'CD1TEST' + '.txt','w')
        
        
    #     def Cd1testReachHeader(xs_number):
                       
    #         adjustable_spacing_1 = ''
            
    #         for i in range(1,5-len(str(self.RivicexsData[xs_number]['Reach']))) : 
                
    #             adjustable_spacing_1 = adjustable_spacing_1 + ' '
            
            
    #         adjustable_spacing_2 = ''
            
    #         for i in range(1,5-len(str(self.RiviceParameters['Rivice XS Interpolation Interval']))) :
                
    #             adjustable_spacing_2 = adjustable_spacing_2 + ' '
            
            
    #         Cd1test.write('PLOT      ELEV                                    INTP')
    #         Cd1test.write('\n')
    #         Cd1test.write('REAC  ' + adjustable_spacing_1 + str(self.RivicexsData[xs_number]['Reach']) + '  REACH FOR WHICH ELEV ARE GIVEN FROM A RELATIVE DATUM')
    #         Cd1test.write('\n')
    #         Cd1test.write('SCAL             1.0       1.0    ' + adjustable_spacing_2 + str(self.RiviceParameters['Rivice XS Interpolation Interval']) + '.')
    #         Cd1test.write('\n')
            
            
    #     def Cd1testXSectionHeader(xs_number):
            
    #         Header_XS_ID = str(self.RivicexsData[xs_number]['XS ID']) + '.'
            
    #         while len(Header_XS_ID) < 10:
                
    #             Header_XS_ID = ' ' + Header_XS_ID
                
                
    #         Header_Distance_XS_precedente = str(self.RivicexsData[xs_number]['Distance XS Precedente']) + '.00'
                
    #         while len(Header_Distance_XS_precedente) < 11 : 
                
    #             Header_Distance_XS_precedente = ' ' + Header_Distance_XS_precedente
                
            
    #         Header_XS = 'SECT' + Header_XS_ID + Header_Distance_XS_precedente
              
    #         Cd1test.write(Header_XS)
    #         Cd1test.write('\n')
            
            
    #     def Cd1testXSectionGeometry(xs_number):
            
            
    #         def length_adjustment(string,desired_length):
                
    #             while len(string) < desired_length:
                    
    #                 string = ' ' + string
                    
    #             return string
    
            
    #         for i in range(0,len(self.RivicexsData[xs_number]['MainChannelGeometry']),6):
                
    #                 ligne = '     '
                
    #                 if i <= len(self.RivicexsData[xs_number]['MainChannelGeometry']) - 2 :
                        
    #                     x1 = str(self.RivicexsData[xs_number]['MainChannelGeometry'][i])
    #                     y1 = str(self.RivicexsData[xs_number]['MainChannelGeometry'][i+1])
                        
    #                     x1 = length_adjustment(x1,10)
    #                     y1 = length_adjustment(y1,10)                        
                        
    #                     ligne = ligne + x1 + y1
                    
    #                 else:
                        
    #                     break
                    
    #                 if i + 2 <= len(self.RivicexsData[xs_number]['MainChannelGeometry']) - 2 :
                        
    #                     x2 = str(self.RivicexsData[xs_number]['MainChannelGeometry'][i+2])
    #                     y2 = str(self.RivicexsData[xs_number]['MainChannelGeometry'][i+3])
                        
    #                     x2 = length_adjustment(x2,10)
    #                     y2 = length_adjustment(y2,10)                        
                        
    #                     ligne = ligne + x2 + y2
                    
    #                 else:
                        
    #                     Cd1test.write(ligne)
    #                     Cd1test.write('\n')
                        
    #                     break
                    
    #                 if i + 4 <= len(self.RivicexsData[xs_number]['MainChannelGeometry']) - 2:
                        
    #                     x3 = str(self.RivicexsData[xs_number]['MainChannelGeometry'][i+4])
    #                     y3 = str(self.RivicexsData[xs_number]['MainChannelGeometry'][i+5])
                        
    #                     x3 = length_adjustment(x3,10)
    #                     y3 = length_adjustment(y3,10)                        
                        
    #                     ligne = ligne + x3 + y3
                        
    #                 else:
                        
    #                     Cd1test.write(ligne)
    #                     Cd1test.write('\n')
                        
    #                     break
                        
    #                 Cd1test.write(ligne)
    #                 Cd1test.write('\n')
            
            

    #     xs_number_precedent = ''
        
    #     for xs_number in self.RivicexsData:
                        
    #         if xs_number == '1':
                
    #             Cd1testReachHeader(xs_number)
    #             Cd1testXSectionHeader(xs_number)
    #             Cd1testXSectionGeometry(xs_number)
                
            
    #         elif self.RivicexsData[xs_number]['Reach'] != self.RivicexsData[xs_number_precedent]['Reach']: 
                
    #             Cd1testReachHeader(xs_number)
    #             Cd1testXSectionHeader(xs_number)
    #             Cd1testXSectionGeometry(xs_number)           
            
            
    #         else:
                
    #             Cd1testXSectionHeader(xs_number)
    #             Cd1testXSectionGeometry(xs_number)
                
                        
    #         xs_number_precedent = xs_number

























def removeOldTIFs(project_path):

    tifs_to_remove = glob.glob(project_path +"\\MonteCarlo\TIF\*.tif")

    for _file in tifs_to_remove:
        os.remove(_file)

    print("Removed all tiffs from previous run if present.\n")














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




# def createTempFile(filepath):
#       directory = os.path.dirname(filepath)
#       filename = os.path.splitext(os.path.basename(filepath))[0]
#       ext = os.path.splitext(os.path.basename(filepath))[1]
#       shutil.copy(filepath, directory+"\\"+filename+"_temp"+ext)
#       return directory+"\\"+filename+"_temp"+ext

# def replaceGeo(filePath,textTofind,textToSub):
#     """
#     Finds a keyword entry in .g0* file and replaces it with supplied value

#     Parameters
#     ----------
#     filePath : TYPE
#     DESCRIPTION.
#     textTofind : TYPE
#     DESCRIPTION.
#     textToSub : TYPE
#     DESCRIPTION.

#     Returns
#     -------
#     None.

#     """

#     file = open(filePath,'r')  # open file handle for read

#     file_temp=createTempFile(filePath)
#     result = open(file_temp, 'w')

#     for line in file:

#         if textTofind in line:
#             line = textToSub+"\n"
#             result.write(line)

#     file.close()
#     result.close()

#     shutil.copyfile(file_temp,filePath)
#     os.remove(file_temp)
#