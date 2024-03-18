
"""
Created by Mathieu Fouquet and Jason Duguay

Date: January 02, 2024
"""

import os
import time
import sys
import numpy as np
import random
import shutil
import subprocess
import pandas as pd
import matplotlib.pyplot as plt

class StochRIVICE():
    
    """
    Converts channel geometry from HECRAS files to format required by RIVICE, creates
	required RIVICE input files, launches necessary RIVICE preprocessing applications 
	(i.e. Cd1xe, Cd2pgm), writes TAPE5.txt, launches RIVICE.

    ...

    """
    
    def __init__(self,stochICE):

        """
        Parameters
        ----------
        stochICE : class instance
            Instance of stochICE with HECRAS reach geometry, flow data, reach parameters
        interInt : integer
            Interpolation spacing between RIVICE cross-sections (in meters)
        """
        
        self.stochICE=stochICE
        
        self.simFolder=self.stochICE.prjDir+'\\RIVICE_simulations'

    
        if not os.path.exists(self.simFolder):
            os.makedirs(self.simFolder)
            print("Created .\RIVICE_simulations")
        


    def make_RIVICE_xs(self): 
        
        print("Creating RIVICE cross-sections ...")
        
        self.xs_data = {}
        
        xs_nb = 1
        xs_prec = ''
        
        #handle bridges
        if self.stochICE.bridge:                        
            
            bridge = 1
            
            for xs in self.stochICE.xs_data:
                
                if xs_prec == '':                
                    self.xs_data[str(xs_nb)] = {}
                    self.xs_data[str(xs_nb)]['Hecras xs'] = xs
                    
                elif bridge <= len(self.stochICE.bridgeData) and self.stochICE.bridgeData[str(bridge)]['chainage'] > float(xs):                    
                    self.xs_data[str(xs_nb)] = {}
                    self.xs_data[str(xs_nb)]['Hecras xs'] = xs_prec
                    
                    xs_nb += 1
                    
                    self.xs_data[str(xs_nb)] = {}
                    self.xs_data[str(xs_nb)]['Hecras xs'] = xs
                    
                    xs_nb += 1
                    
                    self.xs_data[str(xs_nb)] = {}
                    self.xs_data[str(xs_nb)]['Hecras xs'] = xs
                    
                    bridge += 1
                           
                elif self.stochICE.xs_data[xs]['Manning']['val_MAIN'] != self.stochICE.xs_data[xs_prec]['Manning']['val_MAIN'] :                    
                    self.xs_data[str(xs_nb)] = {}
                    self.xs_data[str(xs_nb)]['Hecras xs'] = xs_prec
                    
                    xs_nb += 1
                    
                    self.xs_data[str(xs_nb)] = {}
                    self.xs_data[str(xs_nb)]['Hecras xs'] = xs
                
                else :                    
                    self.xs_data[str(xs_nb)] = {}
                    self.xs_data[str(xs_nb)]['Hecras xs'] = xs
            
                xs_nb += 1        
                xs_prec = xs
                
        else:
            
            for xs in self.stochICE.xs_data:
               
                if xs_prec == '':                
                    self.xs_data[str(xs_nb)] = {}
                    self.xs_data[str(xs_nb)]['Hecras xs'] = xs
                    
                elif self.stochICE.xs_data[xs]['Manning']['val_MAIN'] != self.stochICE.xs_data[xs_prec]['Manning']['val_MAIN'] :                    
                    self.xs_data[str(xs_nb)] = {}
                    self.xs_data[str(xs_nb)]['Hecras xs'] = xs_prec
                    
                    xs_nb += 1
                    
                    self.xs_data[str(xs_nb)] = {}
                    self.xs_data[str(xs_nb)]['Hecras xs'] = xs
                
                else :                    
                    self.xs_data[str(xs_nb)] = {}
                    self.xs_data[str(xs_nb)]['Hecras xs'] = xs
            
                xs_nb += 1        
                xs_prec = xs
            

            
    def get_RIVICE_xs_chainage(self):
        
        print("Getting RIVICE chainages ...")
        
        def round_to_multiple(x, base):
            return base * round(x/base)
        
        flag = False
        reste_arrondi = 0            
        xs_nb_precedent = ''
        xs_nb_precedent_precedent = ''
        
        for xs_nb in self.xs_data:
            
            #p_ stands for precedent
            #pp_ stands for precendent precedent
            
            try:
                #check for empty entries
                p_chainage = self.xs_data[xs_nb_precedent]['Rivice Chainage']
                p_hec_xs = int(self.xs_data[xs_nb_precedent]['Hecras xs'])
                hec_xs = int(self.xs_data[xs_nb]['Hecras xs'])
                pp_chainage = self.xs_data[xs_nb_precedent_precedent]['Rivice Chainage']
                pp_hec_xs = int(self.xs_data[xs_nb_precedent_precedent]['Hecras xs'])
            except KeyError:
                pass

            
            if xs_nb == '1':
                self.xs_data[xs_nb]['Rivice Chainage'] = 0
                
            elif xs_nb == str(len(self.xs_data)):
                self.xs_data[xs_nb]['Rivice Chainage'] = round_to_multiple(p_chainage + (p_hec_xs - hec_xs), self.stochICE.riv_interval)
                 
            elif flag:
                self.xs_data[xs_nb]['Rivice Chainage'] =p_chainage + (p_hec_xs - (reste_arrondi + hec_xs))
                flag = False
                
            elif self.xs_data[xs_nb]['Hecras xs'] == self.xs_data[xs_nb_precedent]['Hecras xs'] :                
                self.xs_data[xs_nb]['Rivice Chainage'] = round_to_multiple(pp_chainage + (pp_hec_xs - p_hec_xs), self.stochICE.riv_interval)
                
                p_chainage = self.xs_data[xs_nb]['Rivice Chainage']   
                
                reste_arrondi = round_to_multiple(pp_chainage + (pp_hec_xs - p_hec_xs), self.stochICE.riv_interval) - (pp_chainage + (pp_hec_xs - p_hec_xs))
                flag = True                                                                                    
            
            else:
                self.xs_data[xs_nb]['Rivice Chainage'] = p_chainage + (p_hec_xs - hec_xs)        
            
            xs_nb_precedent_precedent = xs_nb_precedent            
            xs_nb_precedent = xs_nb
    
     
    def get_RIVICE_xs_manning(self):
        
        print("Getting Manning's coeffients ...")
        
        flag = False
        xs_nb_precedent = ''
        
        for xs_nb in self.xs_data:
            
            xs_manning=self.stochICE.xs_data[self.xs_data[xs_nb]['Hecras xs']]['Manning']['val_MAIN']
            
            if xs_nb == '1':    
                self.xs_data[xs_nb]['Manning'] = xs_manning
            
            elif self.xs_data[xs_nb]['Hecras xs'] == self.xs_data[xs_nb_precedent]['Hecras xs']:
                flag = True
                
            elif flag:
                self.xs_data[xs_nb_precedent]['Manning'] = xs_manning
                self.xs_data[xs_nb]['Manning'] = xs_manning
                flag = False
                
            else:
                self.xs_data[xs_nb]['Manning'] = xs_manning

            xs_nb_precedent = xs_nb
            

            
    def get_RIVICE_xs_geometry(self):
        
        for xs_nb in self.xs_data:            
            xs_geometry=self.stochICE.xs_data[self.xs_data[xs_nb]['Hecras xs']]['MainChannelGeometry']['xy']            
            self.xs_data[xs_nb]['MainChannelGeometry'] = xs_geometry


            
    def assign_RIVICE_Reaches(self):
        
        print("Assigning reaches ...")
        
        Reach = 1
        
        xs_nb_precedent = ''
        
        for xs_nb in self.xs_data:
            
            if xs_nb == '1':            
                self.xs_data[xs_nb]['Reach'] = Reach
                               
            elif self.xs_data[xs_nb]['Manning'] != self.xs_data[xs_nb_precedent]['Manning']:               
                Reach +=1                
                self.xs_data[xs_nb]['Reach'] = Reach
                
            else:                
                self.xs_data[xs_nb]['Reach'] = Reach
                    
            xs_nb_precedent = xs_nb
            

            
    def assign_RIVICE_xs_ID(self):
        
        print("Assigning xs IDs ...")
        
        for xs_nb in self.xs_data:            
            self.xs_data[xs_nb]['XS ID'] = int(xs_nb)*1000
        
        
    def compute_dist_prev_xs(self):
        
        xs_nb_precedent = ''
        
        for xs_nb in self.xs_data:            
            if xs_nb == '1':            
                self.xs_data[xs_nb]['Distance XS Precedente'] = 0            
            else:                
                self.xs_data[xs_nb]['Distance XS Precedente'] = self.xs_data[xs_nb]['Rivice Chainage'] - self.xs_data[xs_nb_precedent]['Rivice Chainage']
                 
            xs_nb_precedent =xs_nb   
                
                
    def get_reach_data(self):
        
        self.reach_data = {}
        
        #Recueil du numero de reach et du chainage de sa limite aval
        
        for i in range(1,len(self.xs_data) + 1):
            
            if i == len(self.xs_data):                
                Reach_number = self.xs_data[str(i)]['Reach']
                Reach_end_chainage = self.xs_data[str(i)]['Rivice Chainage']
                
                self.reach_data[str(Reach_number)] = {}
                self.reach_data[str(Reach_number)]['Reach_number'] = Reach_number
                self.reach_data[str(Reach_number)]['Reach_end_chainage'] = Reach_end_chainage
                               
            else:            
                Reach_number = self.xs_data[str(i)]['Reach']
                Reach_end_chainage = self.xs_data[str(i)]['Rivice Chainage']                
                Next_reach_number = self.xs_data[str(i+1)]['Reach']
                                
                if Reach_number != Next_reach_number:                    
                    self.reach_data[str(Reach_number)] = {}
                    self.reach_data[str(Reach_number)]['Reach_number'] = Reach_number
                    self.reach_data[str(Reach_number)]['Reach_end_chainage'] = Reach_end_chainage

                    
        #Recueil de la longueur de chaque reach et du nombre de XS qu'ils comprennent
        #incluant les XS interpolees            
        for i in range(1,len(self.reach_data)+1):
            
            if i == 1:                
                Reach_length = self.reach_data[str(i)]['Reach_end_chainage']                
                Num_of_XS = int(Reach_length)/self.stochICE.riv_interval + 1
                
                self.reach_data[str(i)]['Reach_length'] = Reach_length                
                self.reach_data[str(i)]['Number_of_XS'] = Num_of_XS
                
            else:                            
                Reach_end_chainage = self.reach_data[str(i)]['Reach_end_chainage']
                Prev_reach_end_chainage = self.reach_data[str(i-1)]['Reach_end_chainage']               
                
                Reach_length = Reach_end_chainage - Prev_reach_end_chainage                
                Num_of_XS = Reach_length/self.stochICE.riv_interval + 1
                
                self.reach_data[str(i)]['Reach_length'] = Reach_length                               
                self.reach_data[str(i)]['Number_of_XS'] = Num_of_XS
                
                
        
    def write_Cd1test(self):    
        
        print("Writing Cd1test.txt ...")
        
        
        Cd1test = open(self.stochICE.prjDir + '/' + 'CD1TEST_intp' + '.txt','w')
        

        
        def write_CD1test_reach_header(xs_nb):
                       
            adjustable_spacing_1 = ''
            
            for i in range(1,5-len(str(self.xs_data[xs_nb]['Reach']))) : 
                
                adjustable_spacing_1 = adjustable_spacing_1 + ' '
            
            
            adjustable_spacing_2 = ''
            
            for i in range(1,5-len(str(self.stochICE.riv_interval))) :
                
                adjustable_spacing_2 = adjustable_spacing_2 + ' '
            
            
            Cd1test.write('PLOT      ELEV                                    INTP')
            Cd1test.write('\n')
            Cd1test.write('REAC  ' + adjustable_spacing_1 + str(self.xs_data[xs_nb]['Reach']) + '  REACH FOR WHICH ELEV ARE GIVEN FROM A RELATIVE DATUM')
            Cd1test.write('\n')
            Cd1test.write('SCAL             1.0       1.0    ' + adjustable_spacing_2 + str(self.stochICE.riv_interval) + '.')
            Cd1test.write('\n')
            

            
        def write_CD1test_xs_header(xs_nb):
            
            Header_XS_ID = str(self.xs_data[xs_nb]['XS ID']) + '.'
            
            while len(Header_XS_ID) < 10:
                
                Header_XS_ID = ' ' + Header_XS_ID
                
                
            Header_Distance_xs_prec = str(self.xs_data[xs_nb]['Distance XS Precedente']) + '.00'
                
            while len(Header_Distance_xs_prec) < 11 : 
                
                Header_Distance_xs_prec = ' ' + Header_Distance_xs_prec
                
            
            Header_XS = 'SECT' + Header_XS_ID + Header_Distance_xs_prec
              
            Cd1test.write(Header_XS)
            Cd1test.write('\n')
            

            
        def write_CD1test_xs_Geometry(xs_nb):
                        
            def length_adjustment(string,desired_length):
                
                while len(string) < desired_length:
                    
                    string = ' ' + string
                    
                return string
                
            for i in range(0,len(self.xs_data[xs_nb]['MainChannelGeometry']),6):
                
                    ligne = '     '
                
                    if i <= len(self.xs_data[xs_nb]['MainChannelGeometry']) - 2 :                        
                        x1 = str(self.xs_data[xs_nb]['MainChannelGeometry'][i])
                        y1 = str(self.xs_data[xs_nb]['MainChannelGeometry'][i+1])
                        
                        x1 = length_adjustment(x1,10)
                        y1 = length_adjustment(y1,10)                        
                        
                        ligne = ligne + x1 + y1
                    
                    else:                        
                        break
                    
                    if i + 2 <= len(self.xs_data[xs_nb]['MainChannelGeometry']) - 2 :                        
                        x2 = str(self.xs_data[xs_nb]['MainChannelGeometry'][i+2])
                        y2 = str(self.xs_data[xs_nb]['MainChannelGeometry'][i+3])
                        
                        x2 = length_adjustment(x2,10)
                        y2 = length_adjustment(y2,10)                        
                        
                        ligne = ligne + x2 + y2
                    
                    else:                        
                        Cd1test.write(ligne)
                        Cd1test.write('\n')                        
                        break
                    
                    if i + 4 <= len(self.xs_data[xs_nb]['MainChannelGeometry']) - 2:                        
                        x3 = str(self.xs_data[xs_nb]['MainChannelGeometry'][i+4])
                        y3 = str(self.xs_data[xs_nb]['MainChannelGeometry'][i+5])
                        
                        x3 = length_adjustment(x3,10)
                        y3 = length_adjustment(y3,10)                        
                        
                        ligne = ligne + x3 + y3
                        
                    else:                        
                        Cd1test.write(ligne)
                        Cd1test.write('\n')                        
                        break
                        
                    Cd1test.write(ligne)
                    Cd1test.write('\n')

        xs_nb_precedent = ''
        
        for xs_nb in self.xs_data:
                        
            if xs_nb == '1':
                
                write_CD1test_reach_header(xs_nb)
                write_CD1test_xs_header(xs_nb)
                write_CD1test_xs_Geometry(xs_nb)
                            
            elif self.xs_data[xs_nb]['Reach'] != self.xs_data[xs_nb_precedent]['Reach']:               
                write_CD1test_reach_header(xs_nb)
                write_CD1test_xs_header(xs_nb)
                write_CD1test_xs_Geometry(xs_nb)           
            
            else:   
                write_CD1test_xs_header(xs_nb)
                write_CD1test_xs_Geometry(xs_nb)
                                       
            xs_nb_precedent = xs_nb
                      
        Cd1test.write("STOP")        
        Cd1test.close()
        
             
    def write_Cd1test_for_DOUT7(self):
        
        """Copies the contents of CD1TEST.txt, but removes the INTP."""
        
        Cd1test = open(self.stochICE.prjDir + '/CD1TEST_intp.txt','r')
        Cd1test_no_INTP = open(self.stochICE.prjDir + '/CD1TEST_no_intp.txt','w')
        
        for line in Cd1test:
            
            if "PLOT" in line:                
                Cd1test_no_INTP.write("PLOT      ELEV                                    ")
                Cd1test_no_INTP.write('\n')
                
            else:                
                Cd1test_no_INTP.write(line)


                 
    def get_water_lvl_TestCd2(self):
        
        self.sim_water_lvl = {}
        
        i = 0
        
        for xs in self.stochICE.xs_data:            
            self.sim_water_lvl[xs] = {}            
            self.sim_water_lvl[xs]['Chainage'] = xs           
            self.sim_water_lvl[xs]['Discharge'] = self.stochICE.open_HECRAS_wse['discharge']            
            self.sim_water_lvl[xs]['Water_lvl_elv'] = self.stochICE.open_HECRAS_wse['wse'][self.stochICE.open_HECRAS_wse['chainage'].index(float(xs))]
            self.sim_water_lvl[xs]['Water_lvl_elv_high'] =  self.sim_water_lvl[xs]['Water_lvl_elv'] + 0.5
            self.sim_water_lvl[xs]['Water_lvl_elv_low'] = self.sim_water_lvl[xs]['Water_lvl_elv'] - 0.5            
            self.sim_water_lvl[xs]['Water_lvl_elv'] = round(self.sim_water_lvl[xs]['Water_lvl_elv'],2)
            self.sim_water_lvl[xs]['Water_lvl_elv_high'] = round(self.sim_water_lvl[xs]['Water_lvl_elv_high'],2)
            self.sim_water_lvl[xs]['Water_lvl_elv_low'] = round(self.sim_water_lvl[xs]['Water_lvl_elv_low'],2)
            
            i += 1
        
    def write_Testcd2(self):
        
        print("Writing Testcd2 ...")        
        Testcd2 = open(self.stochICE.prjDir + '/' + 'TESTCD2' + '.txt','w')


        
        def string_length_adjustment(string,desired_length):            

            while len(string) < desired_length:                
                string = ' ' + string                
            return string
                                
                
        
        def write_Testcd2_reach_header(Reach_number):            

            water_lvl_max = -9999.0
            water_lvl_min = 9999.0
            Manning = 0.0
            Reach_length = -9999.0
            Interp_interval = str(self.stochICE.riv_interval) + '.'
            
            for xs_nb in self.xs_data:
                
                if self.xs_data[xs_nb]['Reach'] == Reach_number:                    
                    if self.sim_water_lvl[self.xs_data[xs_nb]['Hecras xs']]['Water_lvl_elv_high'] > water_lvl_max:
                        water_lvl_max = self.sim_water_lvl[self.xs_data[xs_nb]['Hecras xs']]['Water_lvl_elv_high'] 
                        
                    if self.sim_water_lvl[self.xs_data[xs_nb]['Hecras xs']]['Water_lvl_elv_low'] < water_lvl_min:                        
                        water_lvl_min = self.sim_water_lvl[self.xs_data[xs_nb]['Hecras xs']]['Water_lvl_elv_low']
                        
                    if int(self.xs_data[xs_nb]['Rivice Chainage']) > Reach_length :                         
                        Reach_length = self.xs_data[xs_nb]['Rivice Chainage']                        
                        Manning = self.xs_data[xs_nb]['Manning']
                                               
            water_lvl_max = str(round(water_lvl_max + 0.5,1))
            water_lvl_min = str(round(water_lvl_min - 0.5,1))
            Manning = str(Manning)
            Reach_length = str(Reach_length) + '.'
                                   
            while len(Manning) < 22:
                
                if Manning[0] != '0':                
                    Manning = '0' + Manning
                    
                if len(Manning) == 3:                    
                    Manning = Manning + '00'
                    
                if len(Manning) == 4:                    
                    Manning = Manning + '0'                   
                    
                else:            
                    Manning = Manning + ' '
                        
            Reach_length = string_length_adjustment(Reach_length,13)    
            Interp_interval = string_length_adjustment(Interp_interval,12) 
            water_lvl_max = string_length_adjustment(water_lvl_max,10) 
            water_lvl_min = string_length_adjustment(water_lvl_min,10) 
                                        
            Testcd2.write('REACH ' + str(Reach_number))
            Testcd2.write('\n')
            Testcd2.write('         ' + str(Reach_number) + '    0    0    0    1')
            Testcd2.write('\n')
            Testcd2.write(Manning + '0.0       0.0' + Reach_length + Interp_interval)
            Testcd2.write('\n')
            Testcd2.write('3    ' + water_lvl_min + water_lvl_max)
            Testcd2.write('\n')
            
            
        
        def write_Testcd2_reach_data(Reach_number,Number_of_reaches):
            
            def water_lvl_format_adjusment(string):
                
                if len(string) == 4:                    
                    string = string + '0'                    
                return string
            
            numbering = 1
            
            for xs_nb in self.xs_data:
                
                if self.xs_data[xs_nb]['Reach'] == Reach_number:                                        
                    Reach_xs_nb = str(numbering) 
                    water_lvl_elv = str(self.sim_water_lvl[self.xs_data[xs_nb]['Hecras xs']]['Water_lvl_elv'])   
                    discharge = str(self.sim_water_lvl[self.xs_data[xs_nb]['Hecras xs']]['Discharge'])
                    water_lvl_elv_high = str(self.sim_water_lvl[self.xs_data[xs_nb]['Hecras xs']]['Water_lvl_elv_high'])
                    water_lvl_elv_low = str(self.sim_water_lvl[self.xs_data[xs_nb]['Hecras xs']]['Water_lvl_elv_low'])
                    
                    water_lvl_elv = water_lvl_format_adjusment(water_lvl_elv)
                    water_lvl_elv_high = water_lvl_format_adjusment(water_lvl_elv_high)
                    water_lvl_elv_low = water_lvl_format_adjusment(water_lvl_elv_low)
 
                    Reach_xs_nb = string_length_adjustment(Reach_xs_nb,9)
                    water_lvl_elv = string_length_adjustment(water_lvl_elv,11)
                    discharge = string_length_adjustment(discharge,10)
                    water_lvl_elv_high = string_length_adjustment(water_lvl_elv_high,8)
                    water_lvl_elv_low = string_length_adjustment(water_lvl_elv_low,8)
                    
                    Testcd2.write(Reach_xs_nb + water_lvl_elv + discharge + '                                ' + water_lvl_elv_low + water_lvl_elv_high)
                    Testcd2.write('\n')
             
                numbering += 1
            
            
            if Reach_number < Number_of_reaches:                
                Testcd2.write('\n')
                
            else:               
                Testcd2.write('0')
                
        
        Number_of_reaches = 0
        
        for xs_nb in self.xs_data:            
            if self.xs_data[xs_nb]['Reach'] > Number_of_reaches:                
                Number_of_reaches = self.xs_data[xs_nb]['Reach']
    
    
        for Reach in range(1,Number_of_reaches + 1):            
            write_Testcd2_reach_header(Reach)
            write_Testcd2_reach_data(Reach,Number_of_reaches)
                                       
        Testcd2.close() 
        time.sleep(0.5) #prevents occasional file sharing error              
                    
                
                
    def launch_Cd1xe_no_INTP(self):
        
        print("Running Cd1xe_e.bat ...")
        
        os.popen('copy CD1TEST_no_intp.txt CD1TEST.txt')
        os.startfile(self.stochICE.prjDir + "\\" + "Cd1x_e.bat")  



    def launch_Cd2pgm_a(self):
        
        print("Running Cd2pgm_a.exe ...")
        os.startfile(self.stochICE.prjDir + "\\" + "Cd2pgm_a.exe")     


    def set_time_parameters(self):
        
        self.sim_duration=self.stochICE.riv_days*(60*60*24)
        
        self.nbr_timesteps=int(self.sim_duration/self.stochICE.riv_timestep)
    
    def set_profile_times(self):

        if self.stochICE.riv_days % self.stochICE.riv_profile_interval == 0:
            
            self.profile_times=[]
            
            self.nbr_profiles=int(self.stochICE.riv_days/self.stochICE.riv_profile_interval)
            
            for i in range(self.nbr_profiles):
                
                self.profile_times.append(int((i+1)*self.stochICE.riv_profile_interval*60*60*24)/self.stochICE.riv_timestep)
            
        else:
            print('\nFATAL ERROR: &%"! ... specify a factor of riv_sim_days and then try again ...')
            sys.exit()


    def get_normal_distribution_sample(self,mean,std,samples):
        
        """
        Gets a normal distribution of n samples of the variable centered around its mean, 
        with a standard deviation of std.
    
        Parameters
        ----------
        mean : float
            The mean of the variable's distribution.
        std : float
            The standard deviation of the distribution.
        samples : int
            Number of samples to return.
    
        Returns
        -------
        list
            List of n samples in the distribution.
    
        """
    
        rng = np.random.default_rng()
        values = rng.normal(mean, std, size=samples)
        
        return ["%.2f" % elem for elem in list(values)]

    def init_default_ice_parms(self):
        
        """
        Populates each simulations ice parameters with default values.
        
        Specify the default parameters here. Values changed in
        stochastic analysis overwrite default parameters.
         
        """
        self.sim_data={}

        
        for sim in range(self.stochICE.NSims):
            self.sim_data['sim_%d'%sim]={}
        
            # Ice variables
            self.sim_data['sim_%d'%sim]['DEPOPT']=1
            self.sim_data['sim_%d'%sim]['VDEP']=1.2
            self.sim_data['sim_%d'%sim]['DIAICE']=0.6
            self.sim_data['sim_%d'%sim]['FRMAX']=0.5
            self.sim_data['sim_%d'%sim]['EROPT']=99
            self.sim_data['sim_%d'%sim]['VERODE']=1.5
            self.sim_data['sim_%d'%sim]['FTRLIM']=99
            self.sim_data['sim_%d'%sim]['LEOPT']=3
            self.sim_data['sim_%d'%sim]['Frontthick']=0.1
            self.sim_data['sim_%d'%sim]['VFACTR']=99
            self.sim_data['sim_%d'%sim]['POROSC']=0.5
            self.sim_data['sim_%d'%sim]['POROSFS']=0.5
            self.sim_data['sim_%d'%sim]['SLUSHT']=0.2
            self.sim_data['sim_%d'%sim]['COHESN']=99
            self.sim_data['sim_%d'%sim]['NBRGSW']=1
            self.sim_data['sim_%d'%sim]['RLOCBRG']=1200
            self.sim_data['sim_%d'%sim]['DAYSBR']=0.0
            self.sim_data['sim_%d'%sim]['BRIDTH']=0.3
            self.sim_data['sim_%d'%sim]['THERMD']=0.3
            self.sim_data['sim_%d'%sim]['NSHEDF']=0
            self.sim_data['sim_%d'%sim]['ISTOP']=2
            self.sim_data['sim_%d'%sim]['IPRTYPE']=1
            self.sim_data['sim_%d'%sim]['LIMITOUT']=20
            self.sim_data['sim_%d'%sim]['TIMETRIGINC']=50
            self.sim_data['sim_%d'%sim]['NTOTHER']=1
            self.sim_data['sim_%d'%sim]['SECTIONINC']=100
            self.sim_data['sim_%d'%sim]['NOTSURE1']=499
            self.sim_data['sim_%d'%sim]['NOTSURE2']=483
            self.sim_data['sim_%d'%sim]['ICEGENMETHOD']=2
            self.sim_data['sim_%d'%sim]['HeatLC']=0.0
            self.sim_data['sim_%d'%sim]['IceVol_init']=[(0,0)]
            self.sim_data['sim_%d'%sim]['IceVol_init']=[(2800,0)]
            self.sim_data['sim_%d'%sim]['ZZK1TAN']=0.218
            self.sim_data['sim_%d'%sim]['ZZK2']=7.52
            self.sim_data['sim_%d'%sim]['ICENOPT']=2
            self.sim_data['sim_%d'%sim]['IX']=0 #HAVE TO GET THIS FROM THE SIMULATION
            self.sim_data['sim_%d'%sim]['FACTOR1']=".5"
            self.sim_data['sim_%d'%sim]['FACTOR2']=".12"
            self.sim_data['sim_%d'%sim]['FACTOR3']=".027"
            self.sim_data['sim_%d'%sim]['CNBED']=".026"
            self.sim_data['sim_%d'%sim]['IBORD']=1
            self.sim_data['sim_%d'%sim]['DAYBORDSTART']=31
            self.sim_data['sim_%d'%sim]['BORDUPBRK']=1.0
            self.sim_data['sim_%d'%sim]['BORDWNBRK']=1.0
            self.sim_data['sim_%d'%sim]['BORDCOEF1']=0.1
            self.sim_data['sim_%d'%sim]['BORDCOEF2']=0.15
            self.sim_data['sim_%d'%sim]['BRDTHK']=0.15
            self.sim_data['sim_%d'%sim]['MELTOPT']=1
            self.sim_data['sim_%d'%sim]['MELTSTART']=1000000
            self.sim_data['sim_%d'%sim]['HeatTC']=0.0
            self.sim_data['sim_%d'%sim]['DISTTRIGINC']=50
            self.sim_data['sim_%d'%sim]['NDOTHER']=1
            
            # Boundary conditions
            self.sim_data['sim_%d'%sim]['Q']=0
            self.sim_data['sim_%d'%sim]['DWSE']=0
            self.sim_data['sim_%d'%sim]['IceVol']=0

    def set_stochastic_variables(self):
        
        #Iterate over simulations
        for sim, value in self.sim_data.items():
            
            #Iterate over stochastic variables
            for variable, dist_parms in self.stochICE.riv_stochvars.items():
                
                if variable == 'RLOCBRG':
                    sample=self.get_normal_distribution_sample(dist_parms[0], dist_parms[1], dist_parms[2])
                    self.sim_data[sim][variable]=round(float(random.choice(sample)))
                else:  
                    sample=self.get_normal_distribution_sample(dist_parms[0], dist_parms[1], dist_parms[2])
                    self.sim_data[sim][variable]=float(random.choice(sample))
                
    def make_sim_folders(self):

        #makes simulation directories
        
        
        if not os.path.exists(self.simFolder):
            os.makedirs(self.simFolder)
        
        for sim in range(self.stochICE.NSims):
            
            simPath = self.simFolder+'\\sim_%s' %sim
            
            if not os.path.exists(simPath):
                os.makedirs(simPath)

    def delete_sim_folders(self):
        
        #deletes simulation directories
        for fname in os.listdir(self.simFolder):
            
                folderpath = os.path.join(self.simFolder, fname)
                
                if os.path.isdir(folderpath):
                        if 'sim_' in fname:
                            shutil.rmtree(folderpath)
                            # os.rmdir(folderpath)
    
    def call_write_TAPE5(self):
        
        for sim in range(self.stochICE.NSims):
            self.sim='sim_%d' %sim
            # print(self.sim)
            self.write_TAPE5()

    def write_TAPE5(self):
        
        print("Writing TAPE5 ...")
        
        def string_length_adjustment(string,length,position):
            
            if position == "F":            
                while len(string) < length:                    
                    string = " " + string
            
            elif position == "B":            
                while len(string) < length:                    
                    string = string + " "
                      
            else:               
                print('Error: The position (3rd argument), should be F (for front adjustment) or B (for back adjustment)')                
            return string

        
        """
        Liste de variables locales qui devront être précisées par
        l'utilisateur (pour l'instant ces variables sont "hardcoded")
    
        """
        # Project solution
        Project_title = "RIVICE MODEL"
        
        # Hydraulic solution options (voir le manuel de l'utilisateur de RIVICE pour plus de détails sur ces options)
        IOpt_1 = 1
        IOpt_2 = 2
        IOpt_3 = 1
        IOpt_4 = 1
        IOpt_5 = 2
        KOpt = 2
        
        # Water quality solution options (voir le manuel de l'utilisateur de RIVICE pour plus de détails sur ces options)
        JOpt_1 = 1
        JOpt_2 = 2
        JOpt_3 = 2
        
        # Water quality parameter options (voir le manuel de l'utilisateur de RIVICE pour plus de détails sur ces options)
        NPARM = 1
        
        # Parameter options (voir le manuel de l'utilisateur de RIVICE pour plus de détails sur ces options)
        WQPAR = ["T"]
        INOP = [0]
        OUTOP = [1]
        DOCALC = [0]
        
        # Time Parameters (voir le manuel de l'utilisateur de RIVICE pour plus de détails sur ces options)
        NPER = 1
        
        #--------------------
        """
        JD - The input script will need the user to easily modify the duration and time step.
        Ideally, they will just input the desired duration and time step and NINC will automatically
        be calculated.
        """
        
        NINC = self.nbr_timesteps #(Number of iterations to run, 43200/8640 = 5 second dt)
        PERIOD = self.sim_duration #Simulation will last a half of a day (in seconds)
        #--------------------
        
        RATIO = ""
        MAXITR = "" 
        EPS = ""
        LPER = ""
        NRINC = ""
        
        # Network Parameters (voir le manuel de l'utilisateur de RIVICE pour plus de détails sur ces options)       
        NREACH = self.xs_data[str(len(self.xs_data))]["Reach"] # -> Cette variable n'a pas besoin d'être précisée par l'utilisateur
        NNODE = NREACH + 1 # -> Cette variable n'a pas besoin d'être précisée par l'utilisateur
        NCTR = 0
        
        # Water quality parameter coefficients
        KEY = 1       
        # 0 - Time from the beginning for entry IMC (in hours) 
        # 1 - Ambient temperature (in degrees Fahrenheit) 
        # 2 - Relative humidity (in percent) 
        # 3 - Wind velocity at 2 m (in miles/hour)
        # 4 - Net solar flux (in BTU/ft²/day)
        # 5 - Net atmospheric flux (in BTU/ft²/day) 
        # 6 - Atmospheric pressure (in mm Hg) 
        T_inputs = [[0.,2900.,2910.,9000.],\
                    [32.,32.,32.,32.],\
                    [25.,25.,25.,25.],\
                    [5.,5.,5.,5.],\
                    [0.,0.,0.,0.],\
                    [0.,0.,0.,0.],\
                    [760.,760.,760.,760.,]]
        
        # Upstream discharge boundary condition parameters   
        # Q_bc_inputs = [[0.0],\
        #             [float('{:,.2f}'.format(self.stochICE.flows[0]))]] #min discharge in range   
        Q_bc_inputs = [[0.0],\
                    [float('{:,.2f}'.format(self.sim_data[self.sim]['Q']))]] #min discharge in range     
        Q_bc_Type = 2                        #Boundary condition type
        Q_bc_Time_dep = 1                    #Boundary condition time dependance
        Q_bc_Serie_len = 1 #Boundary condition serie length
        Q_bc_Intr_type = 1                   #Interpolation type used to interpolate between the
                                             #boundary condition serie terms              
                
        # Downstream water surface elevation boundary condition parameters
        # self.ds_wse=float('{:,.2f}'.format(self.stochICE.open_HECRAS_wse['wse'][-1]))
        self.ds_wse=self.sim_data[self.sim]['DWSE']
        H_bc_inputs = [[0.0],\
                    [self.ds_wse]]
        H_bc_Type = 1                        #Boundary condition type
        H_bc_Time_dep = 1                    #Boundary condition time dependance
        H_bc_Serie_len = 1 #Boundary condition serie length
        H_bc_Intr_type = 1                   #Interpolation type used to interpolate between the
                                             #boundary condition serie terms   
                
        # Hydraulic hydrographs and profiles parameters
        """Provides hydrograph at inlet only. If you want more, you must specify their exact chainage. 
        The hdyrograph code in the tape5 write funciton needs work if more hydrographs are required.
        """
        HG_xs_chainage = [0.0] 
        HG_PF_time_step = self.profile_times # Corresponds to the time step number
                                             # where the profiles will be
                                             # printed, not the real time
                                                                                  
        # Upstream water temperature boundary condition parameters
        # 0 - Time step number
        # 1 - Water temperature (Fahrenheit)
        # 2 - Dispersive flux
        # 3 - Total flux       
        T_bc_inputs = [[0.0,1440.0,2880.0],\
                    [32.0,32.0,32.0],\
                    [0.0,0.0,0.0],\
                    [0.0,0.0,0.0]]
        T_bc_Type = 1 #Boundary condition type
        T_bc_Time_dep = 2 #Boundary condition time dependance
        T_bc_Serie_len = len(T_bc_inputs[0]) #Boundary condition serie length
        
        # Water quality graphs and profile output parameters 
        WQG_xs_chainage = [0.0] #***Chainage de XS ou distance de la XS a partir de debut du reach en pieds? Valider cette info
        WQG_PF_time_step = self.profile_times # Corresponds to the time step number
                                               # where the profiles will be
                                               # printed, not the real time
                                               
        # Upstream incoming ice volume boundary condition
        # The last entry must be equal to the number of timesteps (i.e. variable NINC), if not - crash!
        
        """
        Logic for incoming ice volumes
        """
        IV_bc_inputs = []

        if self.stochICE.riv_ice_start != 0:
            
            IV_bc_inputs.append([0,0])
            
            start_of_ice=self.stochICE.riv_ice_start*(60*60*24)
            time_steps=int(start_of_ice/self.stochICE.riv_timestep)
            IV_bc_inputs.append([time_steps,self.sim_data[self.sim]['IceVol']])

        if self.stochICE.riv_ice_start == 0:
            
            IV_bc_inputs.append([0,self.sim_data[self.sim]['IceVol']])

        if self.stochICE.riv_ice_end != self.stochICE.riv_days:
            print("I am in here")
            end_of_ice=self.stochICE.riv_ice_end*(60*60*24)
            time_steps=int(end_of_ice/self.stochICE.riv_timestep)
            IV_bc_inputs.append([time_steps,0])
            
            end_of_ice=self.stochICE.riv_days*(60*60*24)
            time_steps=int(end_of_ice/self.stochICE.riv_timestep)
            IV_bc_inputs.append([time_steps,0])

        if self.stochICE.riv_ice_end == self.stochICE.riv_days:
            
            end_of_ice=self.stochICE.riv_days*(60*60*24)
            time_steps=int(end_of_ice/self.stochICE.riv_timestep)
            IV_bc_inputs.append([time_steps,self.sim_data[self.sim]['IceVol']])

        print(IV_bc_inputs)
                                               
        # RIVICE default parameters (voir le manuel de l'utilisateur de RIVICE pour plus de détails sur ces options)        
        DEPOPT = self.sim_data[self.sim]['DEPOPT']
        VDEP = self.sim_data[self.sim]['VDEP']
        DIAICE = self.sim_data[self.sim]['DIAICE']
        FRMAX = self.sim_data[self.sim]['FRMAX']
        EROPT = self.sim_data[self.sim]['EROPT']
        VERODE = self.sim_data[self.sim]['VERODE']
        FTRLIM = self.sim_data[self.sim]['FTRLIM']
        LEOPT = self.sim_data[self.sim]['LEOPT']
        Frontthick = self.sim_data[self.sim]['Frontthick']
        VFACTR = self.sim_data[self.sim]['VFACTR']
        POROSC = self.sim_data[self.sim]['POROSC']
        POROSFS = self.sim_data[self.sim]['POROSFS']
        SLUSHT = self.sim_data[self.sim]['SLUSHT']
        COHESN = self.sim_data[self.sim]['COHESN']
        NBRGSW = self.sim_data[self.sim]['NBRGSW']
        
        RLOCBRG = self.sim_data[self.sim]['RLOCBRG']
        DAYSBR = self.sim_data[self.sim]['DAYSBR']
        BRIDTH = self.sim_data[self.sim]['BRIDTH']
        THERMD = self.sim_data[self.sim]['THERMD']
        NSHEDF = self.sim_data[self.sim]['NSHEDF']
        ISTOP = self.sim_data[self.sim]['ISTOP'] #Changed this from 0 to 2 to prevent a crash
        IPRTYPE = self.sim_data[self.sim]['IPRTYPE']
        LIMITOUT = self.sim_data[self.sim]['LIMITOUT']
        TIMETRIGINC = self.sim_data[self.sim]['TIMETRIGINC']
        NTOTHER = self.sim_data[self.sim]['NTOTHER']
        SECTIONINC = self.sim_data[self.sim]['SECTIONINC']
        DISTTRIGINC = self.sim_data[self.sim]['DISTTRIGINC']
        NDOTHER = self.sim_data[self.sim]['NDOTHER']
        ICEGENMETHOD = self.sim_data[self.sim]['ICEGENMETHOD']
        Heat_loss_coef = self.sim_data[self.sim]['HeatLC'] #For simple calculation method
        
        # RIVICE INPUT TYPE J parameters (voir le manuel de l'utilisateur de RIVICE pour plus de détails sur ces options)
        ZZK1TAN = self.sim_data[self.sim]['ZZK1TAN']
        ZZK2 = self.sim_data[self.sim]['ZZK2']
        ICENOPT = self.sim_data[self.sim]['ICENOPT']
        FACTOR1 = self.sim_data[self.sim]['FACTOR1']
        FACTOR2 = self.sim_data[self.sim]['FACTOR2']
        FACTOR3 = self.sim_data[self.sim]['FACTOR3']
        CNBED = self.sim_data[self.sim]['CNBED']
        IBORD = self.sim_data[self.sim]['IBORD']
        DAYBORDSTART = self.sim_data[self.sim]['DAYBORDSTART']
        BORDUPBRK = self.sim_data[self.sim]['BORDUPBRK']
        BORDWNBRK = self.sim_data[self.sim]['BORDWNBRK']
        BORDCOEF1 = self.sim_data[self.sim]['BORDCOEF1']
        BORDCOEF2 = self.sim_data[self.sim]['BORDCOEF2']
        BRDTHK = self.sim_data[self.sim]['BRDTHK']
        MELTOPT = self.sim_data[self.sim]['MELTOPT']
        MELTSTART = self.sim_data[self.sim]['MELTSTART']
        Heat_coef_water_to_ice = self.sim_data[self.sim]['HeatTC']

        
        """
        Écriture du fichier de contrôle TAPE5        
        """
        
        TAPE5 = open(self.stochICE.prjDir + '/RIVICE_simulations' + '\\%s' % self.sim + '\\TAPE5.txt','w')
         
        # Writing project name
        if len(Project_title) > 80:            
            print("Error: The project title should be 80 characters or less. Please shorten the project title")
        
        flag = 1
        
        while len(Project_title) < 80:                
            if flag == 1:                
                Project_title = "*" + Project_title               
                flag = 0
                
            else:                
                Project_title = Project_title + "*"                
                flag = 1
                
        TAPE5.write(Project_title)
        TAPE5.write("\n")
        
        # Writing Hydraulic solution options
        Hyd_sol_Options = [IOpt_1,IOpt_2,IOpt_3,IOpt_4,IOpt_5,KOpt]
        
        line = ""
        
        for i in range(len(Hyd_sol_Options)):            
            line = line + "         " + str(Hyd_sol_Options[i])
            
        TAPE5.write(line)
        TAPE5.write("\n")
        
        # Writing Water quality solution options
        WQ_sol_Options = [JOpt_1,JOpt_2,JOpt_3]
        
        line = ""
        
        for i in range(len(WQ_sol_Options)):            
            line = line + "         " + str(WQ_sol_Options[i])
            
        TAPE5.write(line)
        TAPE5.write("\n")
        
        # Writing Water quality parameter options
        line = str(NPARM)
        
        while len(line) < 10:            
            line = " " + line
        
        TAPE5.write(line)
        TAPE5.write("\n")
        
        # Writing Parameter options        
        line = ""
        
        for i in range(len(WQPAR)):            
            while len(WQPAR[i]) < 10:
                
                WQPAR[i] = WQPAR[i] + " "
                
            line = line + WQPAR[i]   
            line = line + "         " + str(INOP[i])
            line = line + "         " + str(OUTOP[i])
            line = line + "         " + str(DOCALC[i])
            
            TAPE5.write(line)
            TAPE5.write("\n")
            
            line = ""
            
        # Writing Time Parameters
        NPER = string_length_adjustment(str(NPER),10,"F")
        NINC = string_length_adjustment(str(NINC),10,"F")
        PERIOD = string_length_adjustment(str(PERIOD),10,"F")
            
        if RATIO != "":        
            print("") # -> Cette partie du script est présentement facultative,
                      #    celle-ci pourra être complétée ultérieurement.
                      
        if MAXITR != "":        
            print("") # -> Cette partie du script est présentement facultative,
                      #    celle-ci pourra être complétée ultérieurement.
                      
        if EPS != "":
            print("") # -> Cette partie du script est présentement facultative,
                      #    celle-ci pourra être complétée ultérieurement.
                      
        if LPER != "":
            print("") # -> Cette partie du script est présentement facultative,
                      #    celle-ci pourra être complétée ultérieurement.
                      
        if NRINC != "":
            print("") # -> Cette partie du script est présentement facultative,
                      #    celle-ci pourra être complétée ultérieurement.
                      
        line = NPER + NINC + PERIOD + RATIO + MAXITR + EPS + LPER + NRINC
        
        TAPE5.write(line)
        TAPE5.write("\n")
        
        # Writing Network Parameters        
        NREACH = string_length_adjustment(str(NREACH),10,"F")
        NNODE = string_length_adjustment(str(NNODE),10,"F")
        NCTR = string_length_adjustment(str(NCTR),10,"F")
        
        line = NREACH + NNODE + NCTR
        
        TAPE5.write(line)
        TAPE5.write("\n")
        
        # Writing  Reach-node Connectivity Data
        # -> Since RIVICE can only model a single channel, reaches informations
        #    are listed for every reach following the reaches numbering order
        for i in range(1,int(NREACH)+1):
            
            line = ""
            
            for j in range(4):                
                if j == 3:                    
                    line = line + string_length_adjustment(str(i + 1),10,"F")
                    
                else:
                    line = line + string_length_adjustment(str(i),10,"F")
                
            TAPE5.write(line)
            TAPE5.write("\n")
            
        # Writing DOUT7 file content       
        Dout7 = open(self.stochICE.prjDir + '/' + 'DOUT7.txt','r')
        
        for line in Dout7:                            
                TAPE5.write(line)
            
        Dout7.close()        
        
        
        # Writing Water quality parameters of the reaches 
        TAPE5.write("Water quality description of the reaches")
        TAPE5.write("\n")
        
        KEY = string_length_adjustment(str(KEY),10,"F")
        
        for i in range(len(WQPAR)):            
            TAPE5.write(WQPAR[i] + KEY)
            TAPE5.write("\n")
            
            if WQPAR[i][0] == "T":            
                T_series = len(T_inputs[0])
                T_series = string_length_adjustment(str(T_series),10,"F")
                TAPE5.write(T_series)
                TAPE5.write("\n")
                          
                for j in range(int(T_series)):                    
                    line = ""
                    
                    for k in range(7):                        
                        param = T_inputs[k][j] 
                        param = string_length_adjustment(str(param),10,"F")
                        
                        line = line + param
                    
                    TAPE5.write(line)
                    TAPE5.write("\n")
                    
            line = ""
              
            for j in range(1,len(self.reach_data)+1):               
                Reach_number = self.reach_data[str(j)]['Reach_number']
                Reach_xs_num = self.reach_data[str(j)]['Number_of_XS']
                
                Reach_number = string_length_adjustment(str(Reach_number),20,'F')
                Reach_xs_num = string_length_adjustment("-" + str(int(Reach_xs_num)),10,'F')
                
                TAPE5.write(Reach_number + Reach_xs_num)
                TAPE5.write("\n")
                
                line = Reach_number + string_length_adjustment("0",10,'F')
                
                TAPE5.write(line)
                TAPE5.write("\n")
                
                TAPE5.write(WQPAR[i] + KEY)
                TAPE5.write("\n")
                
                initial_time = T_inputs[0][0] 
                initial_wat_temp = T_inputs[1][0]
                
                initial_time = string_length_adjustment(str(initial_time),10,"F")
                initial_wat_temp = string_length_adjustment(str(initial_wat_temp),10,"F")
                
                TAPE5.write(initial_time + initial_wat_temp)
                TAPE5.write("\n")
                
                
        # Writing Lateral Inflow parameters of the reaches 
        TAPE5.write("D              LATERAL INFLOW")
        TAPE5.write("\n")
        TAPE5.write("         0")
        TAPE5.write("\n")
        
        # Writing Injection Points parameters of the reaches 
        TAPE5.write("E              INJECTION POINTS")
        TAPE5.write("\n")
        TAPE5.write("         0")
        TAPE5.write("\n")
        
        # Writing Hydraulic Boundary Conditions at each "Node" of the domain 
        TAPE5.write("F              BOUNDARY CONDITIONS")
        TAPE5.write("\n")
        
        reach_num = string_length_adjustment("1",10,'F')
        
        line = reach_num
        
        # Upstream discharge boundary condition definition        
        US_Q_boundary_cdn = [Q_bc_inputs,Q_bc_Type,Q_bc_Time_dep,Q_bc_Serie_len,Q_bc_Intr_type] #US boundary condition parameters
        
        for i in range(1,len(US_Q_boundary_cdn)):
            
            US_Q_boundary_cdn[i] = string_length_adjustment(str(US_Q_boundary_cdn[i]),10,'F')
            
            line = line + US_Q_boundary_cdn[i]
        
        TAPE5.write(line)
        TAPE5.write("\n")
        
        for i in range(len(US_Q_boundary_cdn[0][0])):
            
            time_step = str(US_Q_boundary_cdn[0][0][i])
            Discharge = str(US_Q_boundary_cdn[0][1][i])
            
            time_step = string_length_adjustment(time_step,10,'F')
            Discharge = string_length_adjustment(Discharge,19,'F')

            TAPE5.write(time_step + Discharge)
            TAPE5.write("\n")
        
        # Internal node definition
        for i in range(2,len(self.reach_data)+1):
            
            reach_num = string_length_adjustment(str(i),10,'F')
            
            TAPE5.write(reach_num + "         0")
            TAPE5.write("\n")
        
            
        """
        This is causing a fatal problem. It wrongly makes the reach number equal to 7
        which then crashes RIVICE at line 107 in the TAPE5.txt for the stochICE_RIVICE_testing.py script.
        
        I imagine this logic is useful for multireach channels to specify that the outlet is on the
        nth reach. I think it would be best to have a variable with the number of reaches as an integer
        and then we can just use this interger as reach_num.
        """
        #i += 1 
        #reach_num = string_length_adjustment(str(i),10,'F')
        
        
        reach_num = string_length_adjustment(str('2'),10,'F') #When '1' is used, RIVICE crashes. It only works with '2'. What is this?
        line = reach_num
        
        # Downstream boundary condition definition
        DS_WSE_boundary_cdn = [H_bc_inputs,H_bc_Type,H_bc_Time_dep,H_bc_Serie_len,H_bc_Intr_type]
        
        for i in range(1,len(DS_WSE_boundary_cdn)):
            
            DS_WSE_boundary_cdn[i] = string_length_adjustment(str(DS_WSE_boundary_cdn[i]),10,'F')
            
            line = line + DS_WSE_boundary_cdn[i]
        
        TAPE5.write(line)
        TAPE5.write("\n")
        
        for i in range(len(DS_WSE_boundary_cdn[0][0])):
            
            time_step = str(DS_WSE_boundary_cdn[0][0][i])
            WSE = str(DS_WSE_boundary_cdn[0][1][i])
            
            time_step = string_length_adjustment(time_step,10,'F')
            WSE = string_length_adjustment(WSE,10,'F')

            TAPE5.write(time_step + WSE)
            TAPE5.write("\n")
            
        # Writing Hydraulic Hydrographs and profiles
        TAPE5.write("HYDRAULIC HYDROGRAPHS & PROFILES")
        TAPE5.write("\n")
        
        # Hydrographs parameters
        HG_num = len(HG_xs_chainage)       
        
        TAPE5.write(string_length_adjustment(str(HG_num),10,'F'))
        TAPE5.write("\n")
        
        
        for i in range(HG_num):
            
            for j in range(1,len(self.xs_data)+1):
                
                if int(HG_xs_chainage[i]) == self.xs_data[str(j)]['Rivice Chainage']:
                    
                    reach_num = self.xs_data[str(j)]['Reach']                
                    xs_chainage = HG_xs_chainage[i]
                    
                    reach_num = string_length_adjustment(str(reach_num),10,'F')
                    xs_chainage = string_length_adjustment(str(xs_chainage),10,'F')
                    line = "         1         1         1"                    
                    
                    TAPE5.write(reach_num + xs_chainage + line)
                    TAPE5.write("\n")
                    
                    break
                    
        # Profiles parameters
        PF_num = len(HG_PF_time_step)
        Reach_num = len(self.reach_data)
        
        TAPE5.write(string_length_adjustment(str(PF_num*Reach_num),10,'F'))
        TAPE5.write("\n")
        
        for i in range(PF_num):
            
            time_step = string_length_adjustment(str(int(HG_PF_time_step[i])),10,'F')
            
            for j in range(1,Reach_num+1):
                
                Reach = string_length_adjustment(str(j),10,'F')               
                line = "         1         1         1"
                
                TAPE5.write(Reach + time_step + line)
                TAPE5.write("\n")
                
                
                
        # Writing Water quality boundary conditions
        TAPE5.write("WATER QUALITY BOUNDARY CONDITIONS")
        TAPE5.write("\n")
        
        # Upstream temperature boundary condition definition  
        US_T_boundary_cdn = [T_bc_inputs,T_bc_Type,T_bc_Time_dep,T_bc_Serie_len]
        
        # First node
        line = "         1"
        
        for i in range(1,len(US_T_boundary_cdn)):

            var = US_T_boundary_cdn[i]
            var = string_length_adjustment(str(var),10,'F')
            
            line = line + var
            
        TAPE5.write(line)
        TAPE5.write("\n")
        
        for i in range(len(T_bc_inputs[0])):
            
            line = "T         "
            
            for j in range(len(T_bc_inputs)):
                
                var = T_bc_inputs[j][i]
                var = string_length_adjustment(str(var),10,'F')
                
                line = line + var
                
            TAPE5.write(line)
            TAPE5.write("\n")
            
        # All other nodes
        for i in range(2,len(self.reach_data)+2):
            
            node_num = string_length_adjustment(str(i),10,'F')
            
            TAPE5.write(node_num + "         0")
            TAPE5.write("\n")

    
        # Writing Water quality time graphs and profiles
        TAPE5.write("WATER QUALITY TIME GRAPHS & PROFILES")
        TAPE5.write("\n")

        """
        Note for Matt: This commented section below is not correctly outputting
        the lines in the TAPE5.txt for water quality graphs. It correctly writes
        the number of profiles in the first line, but fails to write the actual
        lines used by RIVICE to output WQ graphs. As I assume these graphs are
        not necessary at the moment, I have hard coded a simple work around 
        below the commented section.
        """

        
        # Time graphs parameters
        # WQG_num = len(WQG_xs_chainage)       
        
        # TAPE5.write(string_length_adjustment(str(WQG_num),10,'F'))
        # TAPE5.write("\n")
        

        # for i in range(WQG_num):
            
        #     for j in range(1,len(self.xs_data)+1):
                
                
        #         if int(WQG_xs_chainage[i]) == self.xs_data[str(j)]['Rivice Chainage']:
                    
        #             reach_num = self.xs_data[str(j)]['Reach']                
        #             xs_chainage = WQG_xs_chainage[i]
                    
        #             reach_num = string_length_adjustment(str(reach_num),10,'F')
        #             xs_chainage = string_length_adjustment(str(xs_chainage),10,'F')
        #             line = "         1"                    
                    
        #             TAPE5.write(reach_num + xs_chainage + line)
        #             TAPE5.write("\n")
                    
        #             break

        """
        This is the hard coded work around which will need to be changed.
        """
        TAPE5.write(string_length_adjustment(str('1'),10,'F')+"\n")
        TAPE5.write("         1" + "       0.0" + "         1" + "\n")

                   
        # Profiles parameters
        PF_num = len(WQG_PF_time_step)
        Reach_num = len(self.reach_data)
        
        TAPE5.write(string_length_adjustment(str(PF_num*Reach_num),10,'F'))
        TAPE5.write("\n")
        
        for i in range(PF_num):
            
            time_step = string_length_adjustment(str(int(WQG_PF_time_step[i])),10,'F')
            
            for j in range(1,Reach_num+1):
                
                Reach = string_length_adjustment(str(j),10,'F')               
                line = "         1"
                
                TAPE5.write(Reach + time_step + line)
                TAPE5.write("\n")

                           
        # Writing Test of RIVICE
        TAPE5.write("******** Test of RIVICE ***********")
        TAPE5.write("\n")
        
        R_def_params = ['DEPOPT',DEPOPT,\
                            'VDEP',VDEP,\
                            'DIAICE',DIAICE,\
                            'FRMAX',FRMAX,\
                            'EROPT',EROPT,\
                            'VERODE',VERODE,\
                            'FTRLIM',FTRLIM,\
                            'LEOPT',LEOPT,\
                            'Frontthick',Frontthick,\
                            'VFACTR',VFACTR,\
                            'POROSC',POROSC,\
                            'POROSFS',POROSFS,\
                            'SLUSHT',SLUSHT,\
                            'COHESN',COHESN,\
                            'NBRGSW',NBRGSW,\
                            'RLOCBRG',RLOCBRG,\
                            'DAYSBR',DAYSBR,\
                            'BRIDTH',BRIDTH,\
                            'THERMD',THERMD,\
                            'NSHEDF',NSHEDF,\
                            'ISTOP',ISTOP,\
                            'IPRTYPE',IPRTYPE,\
                            'LIMITOUT',LIMITOUT,\
                            'TIMETRIGINC',TIMETRIGINC,\
                            'NTOTHER',NTOTHER,\
                            'SECTIONINC',SECTIONINC,\
                            'DISTTRIGINC',DISTTRIGINC,\
                            'NDOTHER',NDOTHER,\
                            'ICEGENMETHOD',ICEGENMETHOD,\
                            'Heat_loss_coef', Heat_loss_coef,\
                            ]
        
        flag_1 = 0
        flag_2 = 0
        flag_3 = 0
        # Rivice default parameters
        for i in range(0,len(R_def_params),2):
            
            name = R_def_params[i]
            value = R_def_params[i+1]
            
            value = string_length_adjustment(str(value),5,'F')
            
            
            if name == '  Frontthick':                                
                
                extra_txt = '  ONLY USED IF LEOPT=3; OTHERWISE IGNORED'
                
                line = value + "   " + name + extra_txt
                
            elif name == 'RLOCBRG':
                
                value_1 = R_def_params[i+1]
                value_2 = R_def_params[i+3]
                value_3 = R_def_params[i+5]
                value_4 = R_def_params[i+7]
                
                value_1 = string_length_adjustment(str(value_1),5,'F')
                value_2 = string_length_adjustment(str(value_2),9,'F')
                value_3 = string_length_adjustment(str(value_3),8,'F')
                value_4 = string_length_adjustment(str(value_4),10,'F')
                
                extra_txt = '             RLOCBRG, DAYSBR, BRIDTH, THERMD'
                
                line = value_1 + value_2 + value_3 + value_4 + extra_txt
                
            elif name == 'DAYSBR' or name == 'BRIDTH' or name == 'THERMD':
                
                flag_3 = 1

            elif name == 'NSHEDF':
                
                extra_txt = '		  Input number of load shed factors'
                
                line = value + "   " + name + extra_txt
                
            elif name == 'ISTOP' :
                
                extra_txt = ' (Cross section number at which program stops)'
                
                line = value + "               " + name + extra_txt
                
            elif name == 'IPRTYPE':
            
                value_1 = R_def_params[i+1]
                value_2 = R_def_params[i+3]
                
                value_1 = string_length_adjustment(str(value_1),5,'F')
                value_2 = string_length_adjustment(str(value_2),5,'F')
                
                extra_txt = "IPRTYPE: 0(GEN'L ONLY),1 (GEN'L + DET'D CALCS);  LIMITOUT"
                
                line = value_1 + value_2 + "          " + extra_txt
                
            elif name == 'LIMITOUT':
                
                flag_3 = 1            
                
            elif name == 'TIMETRIGINC':
                
                value_1 = R_def_params[i+1]
                value_2 = R_def_params[i+3]
                
                value_1 = string_length_adjustment(str(value_1),5,'F')
                value_2 = string_length_adjustment(str(value_2),5,'F')
                
                extra_txt = 'TIMETRIGINC, NTOTHER'
                
                line = value_1 + value_2 + "          " + extra_txt
                
            elif name == 'NTOTHER':
                
                flag_3 = 1 
                
            elif name == 'SECTIONINC':
                
                extra_txt = '               SECTIONINC'
                
                line = value + extra_txt

                flag_1 = 1   
                
            elif name == 'DISTTRIGINC':
                
                value_1 = R_def_params[i+1]
                value_2 = R_def_params[i+3]
                
                value_1 = string_length_adjustment(str(value_1),5,'F')
                value_2 = string_length_adjustment(str(value_2),5,'F')
                
                extra_txt = 'TIMETRIGINC, NTOTHER'
                
                line = value_1 + value_2 + "          " + extra_txt
                
                flag_2 = 1
                
            elif name == 'NDOTHER':
                
                flag_3 = 1 
                
            elif name == 'ICEGENMETHOD':
                
                extra_txt = 'ICEGENMETHOD   (1-DETAILED;2-SIMPLIFIED)'
                
                line = value + '                      ' + extra_txt
            
            elif name == 'Heat_loss_coef':
                
                value_1 = R_def_params[i+1] 
                
                value_1 = string_length_adjustment(str(value_1),4,'F')
                
                extra_txt = 'HEAT LOSS COEFFIFICIENT FOR SIMPLE CALC METHOD'
                
                line = value_1 + "                       " + extra_txt
                
            else:
            
                line = value + "   " + name
                
           
            if flag_1 == 1:
                
                TAPE5.write(line)
                TAPE5.write("\n")
                
                TAPE5.write('     499')
                TAPE5.write("\n")
                
                flag_1 = 0
                
            elif flag_2 == 1:
                
                TAPE5.write(line)
                TAPE5.write("\n")
                
                TAPE5.write('     483')
                TAPE5.write("\n")
                
                flag_2 = 0
                
            elif flag_3 == 1:

                flag_3 = 0
                  
            else:
                TAPE5.write(line)
                TAPE5.write("\n")
            
            
        # Upstream incoming ice volume  
        for i in range(len(IV_bc_inputs)):
            print(IV_bc_inputs[i])
            
            time_step = int(IV_bc_inputs[i][0])
            ice_volume = IV_bc_inputs[i][1]
            
            time_step = string_length_adjustment(str(time_step),10,'F')
            ice_volume = string_length_adjustment(str(ice_volume),8,'F')
            
            extra_txt = '         INCOMING ICE VOLUME FOR EACH TIME STEP'
            
            line = time_step + ice_volume + extra_txt
            
            TAPE5.write(line)
            TAPE5.write("\n")
            
        
        # Writing ICE INFORMATION THRU RIVINH  INPUT TYPE J
        TAPE5.write('ICE INFORMATION THRU RIVINH  INPUT TYPE J *********************************')
        TAPE5.write("\n")
        
        nb_tot_XS = 0
        
        for i in range(1,len(self.reach_data)+1):
            
            reach_XS_nb = int(self.reach_data[str(i)]['Number_of_XS'])
            
            nb_tot_XS = nb_tot_XS + reach_XS_nb
        
        
        value_1 = string_length_adjustment(str(ZZK1TAN),5,'F')
        value_2 = string_length_adjustment(str(ZZK2),10,'F')
        extra_txt = '          ZZK1TAN,  ZZK2            '
        line = value_1 + value_2 + extra_txt
        TAPE5.write(line)
        TAPE5.write("\n")
        
        value_1 = string_length_adjustment(str(ICENOPT),5,'F')
        extra_txt = ' ICENOPT(1-BELTAOS, 2-KGS, 3- USER-DEFINED)   '
        line = value_1 + extra_txt
        TAPE5.write(line)
        TAPE5.write("\n")
        
        value_1 = string_length_adjustment(str(nb_tot_XS),5,'F')
        value_2 = string_length_adjustment(FACTOR1,5,'F')
        value_3 = string_length_adjustment(FACTOR2,5,'F')
        value_4 = string_length_adjustment(str(FACTOR3),5,'F')
        value_5 = string_length_adjustment(str(CNBED),5,'F')
        extra_txt = '          IX,FACTOR1, FACTOR2, FACTOR3, CNBED'
        line = value_1 + value_2 + value_3 + value_4 + value_5 + extra_txt
        print(FACTOR1,FACTOR2)
        TAPE5.write(line)
        TAPE5.write("\n")
        
        value_1 = string_length_adjustment(str(IBORD),5,'F')
        extra_txt = '                 IBORD (1-USER, 2-NEWBURY, 3- MATOUSEK)       '
        line = value_1 + extra_txt
        TAPE5.write(line)
        TAPE5.write("\n")
        
        value_1 = string_length_adjustment(str(DAYBORDSTART),4,'F')
        extra_txt = '                                DAYBORDSTART'
        line = value_1 + extra_txt
        TAPE5.write(line)
        TAPE5.write("\n")
        
        value_1 = string_length_adjustment(str(BORDUPBRK),5,'F')
        value_2 = string_length_adjustment(str(BORDWNBRK),5,'F')
        extra_txt = '                 BORDUPBRK   BORDWNBRK           '
        line = value_1 + value_2 + extra_txt
        TAPE5.write(line)
        TAPE5.write("\n")
        
        value_1 = string_length_adjustment(str(nb_tot_XS),5,'F')
        value_2 = string_length_adjustment(str(BORDCOEF1),5,'F')
        value_3 = string_length_adjustment(str(BORDCOEF2),5,'F')
        value_4 = string_length_adjustment(str(BRDTHK),5,'F')
        extra_txt = '  IX  BORDCOEF1  BORDCOEF2  BRDTHK'
        line = value_1 + value_2 + value_3 + value_4 + extra_txt
        TAPE5.write(line)
        TAPE5.write("\n")
        
        value_1 = string_length_adjustment(str(MELTOPT),5,'F')
        value_2 = string_length_adjustment(str(MELTSTART),22,'B')
        extra_txt = 'MELTOPT (1- uSER COEFF,2- RIVICE ALGO.);MELTSTART'
        line = value_1 + " " + value_2 + extra_txt
        TAPE5.write(line)
        TAPE5.write("\n")
        
        value_1 = string_length_adjustment(str(Heat_coef_water_to_ice),6,'F')
        extra_txt = '                      HEAT TRANSFER COEFFICIENT WATER TO ICE(btu/M2/DAY)'
        line = value_1 + extra_txt
        TAPE5.write(line)
        TAPE5.write("\n")
        TAPE5.close()
        
        
    def launch_Cd1xe_with_INTP(self):
        
        os.popen('copy CD1TEST_intp.txt CD1TEST.txt')
        time.sleep(0.5) #unfortunately necessary to be able to open CD1TEST.txt in next line
        os.startfile(self.stochICE.prjDir + "\\" + "Cd1x_e.bat")
        time.sleep(0.5)
            
        
    def launch_RIVICE(self):
        
        os.startfile(self.stochICE.prjDir + "\\" + "rivice.bat")        
        
       
        
    def clean_RIVICE_case(self):   


        files=['CD1TEST.txt',
               'CD1TEST_intp.txt',
               'CD1TEST_no_intp.txt',
               'DOUT1.TXT',
               'DOUT2.TXT',
               'DOUT7.TXT',
               'INPT10.TXT',
               'INPT11.TXT',
               'INPT101.TXT',
               'INPT111.TXT',
               'psplot.ps',
               'TAPE5.txt', 'TAPE6.txt','TAPE10.txt',
               'TAPE12.txt','TAPE13.txt','TAPE14.txt',
               'TAPE15.txt','TAPE16.txt','TAPE17.txt',
               'TAPE66.txt','TAPE70.txt','TAPE71.txt',
               'TAPE72.txt','TAPE73.txt','TAPE5.txt',
               'intpxs.txt','intpxs1.txt',
               'TESTCD2.txt']

        for file in files:
        
            try:
                os.remove(self.stochICE.prjDir +'\\' + file)
                print('Removed %s file' %file)
            except FileNotFoundError:
                pass
                
        
    def launch_RIVICE_parallel(self):
    
        if self.stochICE.NSims % self.stochICE.NProcs == 0:
            
            group_size=int(self.stochICE.NSims/self.stochICE.NProcs)
            sim_list=list(self.sim_data.keys())
            sim_groups=[sim_list[i:i + group_size] for i in range(0, len(sim_list), group_size)]
            
        else:
            print("Please make NSims/NProcs a whole number.")
            pass
        
        for number, group in enumerate(sim_groups):
            
            processes=[]
            
            print('---- Running RIVICE in parallel on group %s ----' %str(number+1))
            
            for sim in group:
                
                shutil.copy('rivice.bat',self.stochICE.prjDir+'\RIVICE_simulations\%s'%sim)
                shutil.copy('Rivice_Aug6_11d.exe',self.stochICE.prjDir+'\RIVICE_simulations\%s'%sim)
                shutil.copy('Rivice_Aug6_11e.exe',self.stochICE.prjDir+'\RIVICE_simulations\%s'%sim)
                shutil.copy('lf90.eer',self.stochICE.prjDir+'\RIVICE_simulations\%s'%sim)
                shutil.copy('intpxs.txt',self.stochICE.prjDir+'\RIVICE_simulations\%s'%sim)
                shutil.copy('intpxs1.txt',self.stochICE.prjDir+'\RIVICE_simulations\%s'%sim)        
        
            for sim in group:
                
                time.sleep(0.5)
                P=subprocess.Popen([self.stochICE.prjDir+'\RIVICE_simulations\%s' %sim + "\\" + "rivice.bat"], cwd=self.stochICE.prjDir+'\RIVICE_simulations\%s' %sim)
                processes.append(P)
                
            for p in processes:
                
                poll=p.poll()
                while poll is None:
                    time.sleep(1)
                    poll=p.poll()
                    print('Group %s is processing ...' %str(number+1))
            
            # print('Removing unnecessary files from simulation folders.')
            # for sim in group:
                
            #     files=['TAPE6.txt','TAPE10.txt',
            #             'TAPE12.txt','TAPE13.txt','TAPE14.txt',
            #             'TAPE15.txt','TAPE16.txt','TAPE17.txt',
            #             'TAPE66.txt','TAPE70.txt','TAPE71.txt',
            #             'TAPE72.txt','TAPE73.txt',
            #             'intpxs.txt','intpxs1.txt',
            #             'Rivice_Aug6_11d.exe','Rivice_Aug6_11e.exe',
            #             'rivice.bat','lf90.eer',]
                
                
            #     for file in files:
                    
            #         try:
            #             os.remove(self.stochICE.prjDir+'\RIVICE_simulations\%s'%sim +'\\' + file)
            #         except FileNotFoundError:
            #             pass        
        

    def get_profiles(self):
        
        last_profile="P_01_%s" %str(round(self.profile_times[-1])).zfill(6) + ".txt"
        
        self.sim_profiles={}
        
        for sim, value in self.sim_data.items():
            
            try:
            
                self.sim_profiles[sim]=pd.read_csv(self.stochICE.prjDir+'\RIVICE_simulations\%s\%s'%(sim,last_profile), header = None, delimiter= '\s+', index_col=False) 
                self.sim_profiles[sim].columns =['chainage','wse','depth','discharge','velocity','n','C','area','rH','thick']
                max_chainage=self.sim_profiles[sim]['chainage'].max()
                self.sim_profiles[sim]['RIVICE_chainage']=max_chainage-self.sim_profiles[sim]['chainage']
            except FileNotFoundError:
                continue
        
    def plot_profiles(self):
    

        self.get_profiles()
        
        fig, ax = plt.subplots(1,1,frameon=False,constrained_layout=True)
        fig.set_figwidth(6,forward=True)
        fig.set_figheight(3,forward=True)
        
        for sim, profiles in self.sim_profiles.items():
            ax.plot(max(profiles['chainage']) - profiles['chainage'],profiles['wse'],c='darkblue')
            
        ax.plot(max(profiles['chainage']) - profiles['chainage'],profiles['wse']-profiles['depth'],c='darkgray')
        ax.set_xlim([min(profiles['chainage']),max(profiles['chainage'])])   
        
        ax.set_xlabel('chainage (m)', fontsize=10)
        ax.set_ylabel('water surface elevation (m)', fontsize=10)        
        
        
    def plot_prob_exceedance(self,chainage):
        
    
        fig, ax = plt.subplots(1,1,frameon=False,constrained_layout=True)
        fig.set_figwidth(6,forward=True)
        fig.set_figheight(3,forward=True)
        
        water_surface_levels=[]
        probabilities=[]
        
        #get water surface levels at x-section closest to specified chainage
        for sim, notused in self.sim_profiles.items():
            
            profile=self.sim_profiles[sim]
            row=profile['RIVICE_chainage'].sub(chainage).abs().idxmin()
            closest_chainage=profile.iloc[row]['RIVICE_chainage']
            wse=profile.iloc[row]['wse']
            water_surface_levels.append(wse)
            
        #sort low to high
        water_surface_levels.sort()
        
        #rank as probability of exceedance
        for count,wse in enumerate(water_surface_levels):
            probabilities.append(1-(count+1)/(len(water_surface_levels)))
            
        ax.scatter(probabilities,water_surface_levels,c='black',marker='o',s=0.5)
        
        plt.gca().invert_xaxis()
        ax.set_xscale('log')
        fig.suptitle('Chainage: %d m'%closest_chainage)
        ax.set_xlabel('exceedence probability (m)', fontsize=10)
        ax.set_ylabel('water surface elevation (m)', fontsize=10)        
            
        
        
        
        
        
        
            
            
            
            
            
            
            
            
            
            
            
            