# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 13:52:46 2023

@author: dugj2403
"""

import os

class StochRIVICE():
    
    # pass-in stochICE instance 
    def __init__(self,stochICE,interInt):
        
        self.stochICE=stochICE
        self.interInt=interInt

           
    def create_RIVICE_xs(self): 
        
        self.riv_xs_data = {}
        
        xs_number = 1
        xs_prec = ''
        
        if self.stochICE.bridge:                        
            
            bridge_number = 1
            
            for xs in self.stochICE.xs_data:
                
                if xs_prec == '':
                
                    self.riv_xs_data[str(xs_number)] = {}
                    self.riv_xs_data[str(xs_number)]['Hecras xs'] = xs
                    
                    
                elif bridge_number <= len(self.stochICE.bridgeData) and self.stochICE.bridgeData[str(bridge_number)]['chainage'] > float(xs):
                    
                    self.riv_xs_data[str(xs_number)] = {}
                    self.riv_xs_data[str(xs_number)]['Hecras xs'] = xs_prec
                    
                    xs_number += 1
                    
                    self.riv_xs_data[str(xs_number)] = {}
                    self.riv_xs_data[str(xs_number)]['Hecras xs'] = xs
                    
                    xs_number += 1
                    
                    self.riv_xs_data[str(xs_number)] = {}
                    self.riv_xs_data[str(xs_number)]['Hecras xs'] = xs
                    
                    bridge_number += 1
                           
                    
                elif self.stochICE.xs_data[xs]['Manning']['val_MAIN'] != self.stochICE.xs_data[xs_prec]['Manning']['val_MAIN'] :
                    
                    self.riv_xs_data[str(xs_number)] = {}
                    self.riv_xs_data[str(xs_number)]['Hecras xs'] = xs_prec
                    
                    xs_number += 1
                    
                    self.riv_xs_data[str(xs_number)] = {}
                    self.riv_xs_data[str(xs_number)]['Hecras xs'] = xs
                
                else :
                    
                    self.riv_xs_data[str(xs_number)] = {}
                    self.riv_xs_data[str(xs_number)]['Hecras xs'] = xs
            
                xs_number += 1        
                xs_prec = xs
                
        else:
            
            for xs in self.stochICE.xs_data:
                
                if xs_prec == '':
                
                    self.riv_xs_data[str(xs_number)] = {}
                    self.riv_xs_data[str(xs_number)]['Hecras xs'] = xs
                    
                           
                elif self.stochICE.xs_data[xs]['Manning']['val_MAIN'] != self.stochICE.xs_data[xs_prec]['Manning']['val_MAIN'] :
                    
                    self.riv_xs_data[str(xs_number)] = {}
                    self.riv_xs_data[str(xs_number)]['Hecras xs'] = xs_prec
                    
                    xs_number += 1
                    
                    self.riv_xs_data[str(xs_number)] = {}
                    self.riv_xs_data[str(xs_number)]['Hecras xs'] = xs
                
                else :
                    
                    self.riv_xs_data[str(xs_number)] = {}
                    self.riv_xs_data[str(xs_number)]['Hecras xs'] = xs
            
                xs_number += 1        
                xs_prec = xs
            
            
    def compute_RIVICE_xs_chainage(self):
        
        def round_to_multiple(x, base):
            return base * round(x/base)
        
        flag = False
        reste_arrondi = 0            
        xs_number_precedent = ''
        xs_number_precedent_precedent = ''
        
        for xs_number in self.riv_xs_data:
            
            if xs_number == '1':
                self.riv_xs_data[xs_number]['Rivice Chainage'] = 0
                
            elif xs_number == str(len(self.riv_xs_data)):
                self.riv_xs_data[xs_number]['Rivice Chainage'] = round_to_multiple(self.riv_xs_data[xs_number_precedent]['Rivice Chainage'] + (int(self.riv_xs_data[xs_number_precedent]['Hecras xs']) - int(self.riv_xs_data[xs_number]['Hecras xs'])), self.interInt)
                 
            elif flag:
                self.riv_xs_data[xs_number]['Rivice Chainage'] = self.riv_xs_data[xs_number_precedent]['Rivice Chainage'] + (int(self.riv_xs_data[xs_number_precedent]['Hecras xs']) - (reste_arrondi + int(self.riv_xs_data[xs_number]['Hecras xs'])))
                flag = False
                
            elif self.riv_xs_data[xs_number]['Hecras xs'] == self.riv_xs_data[xs_number_precedent]['Hecras xs'] :
                
                self.riv_xs_data[xs_number]['Rivice Chainage'] = round_to_multiple(self.riv_xs_data[xs_number_precedent_precedent]['Rivice Chainage'] + (int(self.riv_xs_data[xs_number_precedent_precedent]['Hecras xs']) - int(self.riv_xs_data[xs_number_precedent]['Hecras xs'])), self.interInt)
                
                self.riv_xs_data[xs_number_precedent]['Rivice Chainage'] = self.riv_xs_data[xs_number]['Rivice Chainage']   
                
                reste_arrondi = round_to_multiple(self.riv_xs_data[xs_number_precedent_precedent]['Rivice Chainage'] + (int(self.riv_xs_data[xs_number_precedent_precedent]['Hecras xs']) - int(self.riv_xs_data[xs_number_precedent]['Hecras xs'])), self.interInt) - (self.riv_xs_data[xs_number_precedent_precedent]['Rivice Chainage'] + (int(self.riv_xs_data[xs_number_precedent_precedent]['Hecras xs']) - int(self.riv_xs_data[xs_number_precedent]['Hecras xs'])))
                flag = True                                                                                    
            
            else:
                self.riv_xs_data[xs_number]['Rivice Chainage'] = self.riv_xs_data[xs_number_precedent]['Rivice Chainage'] + (int(self.riv_xs_data[xs_number_precedent]['Hecras xs']) - int(self.riv_xs_data[xs_number]['Hecras xs']))        
            
            xs_number_precedent_precedent = xs_number_precedent
            
            xs_number_precedent = xs_number
    
     
    def get_RIVICE_xs_manning(self):
        
        flag = False
        xs_number_precedent = ''
        
        for xs_number in self.riv_xs_data:
            
            if xs_number == '1':    
                self.riv_xs_data[xs_number]['Manning'] = self.stochICE.xs_data[self.riv_xs_data[xs_number]['Hecras xs']]['Manning']['val_MAIN']
            
            elif self.riv_xs_data[xs_number]['Hecras xs'] == self.riv_xs_data[xs_number_precedent]['Hecras xs']:
                flag = True
                
            elif flag:
                self.riv_xs_data[xs_number_precedent]['Manning'] = self.stochICE.xs_data[self.riv_xs_data[xs_number]['Hecras xs']]['Manning']['val_MAIN']
                self.riv_xs_data[xs_number]['Manning'] = self.stochICE.xs_data[self.riv_xs_data[xs_number]['Hecras xs']]['Manning']['val_MAIN']
                flag = False
                
            else:
                self.riv_xs_data[xs_number]['Manning'] = self.stochICE.xs_data[self.riv_xs_data[xs_number]['Hecras xs']]['Manning']['val_MAIN']

            xs_number_precedent = xs_number
            
            
    def get_RIVICE_xs_geometry(self):
        
        for xs_number in self.riv_xs_data:
            self.riv_xs_data[xs_number]['MainChannelGeometry'] = self.stochICE.xs_data[self.riv_xs_data[xs_number]['Hecras xs']]['MainChannelGeometry']['xy']
            
            
    def assign_RIVICE_Reaches(self):
        
        Reach = 1
        
        xs_number_precedent = ''
        
        for xs_number in self.riv_xs_data:
            
            if xs_number == '1':
            
                self.riv_xs_data[xs_number]['Reach'] = Reach
                
                
            elif self.riv_xs_data[xs_number]['Manning'] != self.riv_xs_data[xs_number_precedent]['Manning']:
                
                Reach +=1
                
                self.riv_xs_data[xs_number]['Reach'] = Reach
                
            else:
                
                self.riv_xs_data[xs_number]['Reach'] = Reach
        
            
            xs_number_precedent = xs_number
            
            
    def assign_RIVICE_xs_ID(self):
        
        for xs_number in self.riv_xs_data:
            
            self.riv_xs_data[xs_number]['XS ID'] = int(xs_number)*1000
        
        
        
        
    def compute_dist_prev_xs(self):
        
        xs_number_precedent = ''
        
        for xs_number in self.riv_xs_data:
            
            if xs_number == '1':
            
                self.riv_xs_data[xs_number]['Distance XS Precedente'] = 0
            
            else:
                
                self.riv_xs_data[xs_number]['Distance XS Precedente'] = self.riv_xs_data[xs_number]['Rivice Chainage'] - self.riv_xs_data[xs_number_precedent]['Rivice Chainage']
                
                
            xs_number_precedent =xs_number   
                
                
   

    def write_Cd1test(self):    
        
        
        Cd1test = open(self.stochICE.prjDir + '/' + 'CD1TEST' + '.txt','w')
        
        
        def write_CD1test_reach_header(xs_number):
                       
            adjustable_spacing_1 = ''
            
            for i in range(1,5-len(str(self.riv_xs_data[xs_number]['Reach']))) : 
                
                adjustable_spacing_1 = adjustable_spacing_1 + ' '
            
            
            adjustable_spacing_2 = ''
            
            for i in range(1,5-len(str(self.interInt))) :
                
                adjustable_spacing_2 = adjustable_spacing_2 + ' '
            
            
            Cd1test.write('PLOT      ELEV                                    INTP')
            Cd1test.write('\n')
            Cd1test.write('REAC  ' + adjustable_spacing_1 + str(self.riv_xs_data[xs_number]['Reach']) + '  REACH FOR WHICH ELEV ARE GIVEN FROM A RELATIVE DATUM')
            Cd1test.write('\n')
            Cd1test.write('SCAL             1.0       1.0    ' + adjustable_spacing_2 + str(self.interInt) + '.')
            Cd1test.write('\n')
            
            
        def write_CD1test_xs_header(xs_number):
            
            Header_XS_ID = str(self.riv_xs_data[xs_number]['XS ID']) + '.'
            
            while len(Header_XS_ID) < 10:
                
                Header_XS_ID = ' ' + Header_XS_ID
                
                
            Header_Distance_xs_prec = str(self.riv_xs_data[xs_number]['Distance XS Precedente']) + '.00'
                
            while len(Header_Distance_xs_prec) < 11 : 
                
                Header_Distance_xs_prec = ' ' + Header_Distance_xs_prec
                
            
            Header_XS = 'SECT' + Header_XS_ID + Header_Distance_xs_prec
              
            Cd1test.write(Header_XS)
            Cd1test.write('\n')
            
            
        def write_CD1test_xs_Geometry(xs_number):
            
            
            def length_adjustment(string,desired_length):
                
                while len(string) < desired_length:
                    
                    string = ' ' + string
                    
                return string
    
            
            for i in range(0,len(self.riv_xs_data[xs_number]['MainChannelGeometry']),6):
                
                    ligne = '     '
                
                    if i <= len(self.riv_xs_data[xs_number]['MainChannelGeometry']) - 2 :
                        
                        x1 = str(self.riv_xs_data[xs_number]['MainChannelGeometry'][i])
                        y1 = str(self.riv_xs_data[xs_number]['MainChannelGeometry'][i+1])
                        
                        x1 = length_adjustment(x1,10)
                        y1 = length_adjustment(y1,10)                        
                        
                        ligne = ligne + x1 + y1
                    
                    else:
                        
                        break
                    
                    if i + 2 <= len(self.riv_xs_data[xs_number]['MainChannelGeometry']) - 2 :
                        
                        x2 = str(self.riv_xs_data[xs_number]['MainChannelGeometry'][i+2])
                        y2 = str(self.riv_xs_data[xs_number]['MainChannelGeometry'][i+3])
                        
                        x2 = length_adjustment(x2,10)
                        y2 = length_adjustment(y2,10)                        
                        
                        ligne = ligne + x2 + y2
                    
                    else:
                        
                        Cd1test.write(ligne)
                        Cd1test.write('\n')
                        
                        break
                    
                    if i + 4 <= len(self.riv_xs_data[xs_number]['MainChannelGeometry']) - 2:
                        
                        x3 = str(self.riv_xs_data[xs_number]['MainChannelGeometry'][i+4])
                        y3 = str(self.riv_xs_data[xs_number]['MainChannelGeometry'][i+5])
                        
                        x3 = length_adjustment(x3,10)
                        y3 = length_adjustment(y3,10)                        
                        
                        ligne = ligne + x3 + y3
                        
                    else:
                        
                        Cd1test.write(ligne)
                        Cd1test.write('\n')
                        
                        break
                        
                    Cd1test.write(ligne)
                    Cd1test.write('\n')
            
            

        xs_number_precedent = ''
        
        for xs_number in self.riv_xs_data:
                        
            if xs_number == '1':
                
                write_CD1test_reach_header(xs_number)
                write_CD1test_xs_header(xs_number)
                write_CD1test_xs_Geometry(xs_number)
                
            
            elif self.riv_xs_data[xs_number]['Reach'] != self.riv_xs_data[xs_number_precedent]['Reach']: 
                
                write_CD1test_reach_header(xs_number)
                write_CD1test_xs_header(xs_number)
                write_CD1test_xs_Geometry(xs_number)           
            
            
            else:
                
                write_CD1test_xs_header(xs_number)
                write_CD1test_xs_Geometry(xs_number)
                
                        
            xs_number_precedent = xs_number
            
            
        Cd1test.write("STOP")
        
        Cd1test.close()
        
        
        
    def write_Cd1test_for_DOUT7(self):
        
        if not os.path.exists(self.stochICE.prjDir + '/DOUT7'):
            
            os.makedirs(self.stochICE.prjDir + '/DOUT7')
            
            
        Cd1test = open(self.stochICE.prjDir + '/CD1TEST.txt','r')
        
        Cd1test_no_INTP = open(self.stochICE.prjDir + '/DOUT7/CD1TEST.txt','w')
        
        
        for line in Cd1test:
            
            if "PLOT" in line:
                
                Cd1test_no_INTP.write("PLOT      ELEV                                    ")
                Cd1test_no_INTP.write('\n')
                
            else:
                
                Cd1test_no_INTP.write(line)
                
        
            
    def get_water_lvl_TestCd2(self):
        
        self.sim_water_lvl = {}
        
        for xs in self.stochICE.xs_data:
            
            self.sim_water_lvl[xs] = {}
            self.sim_water_lvl[xs]['Chainage'] = xs
            self.sim_water_lvl[xs]['Discharge'] = 125.0
            self.sim_water_lvl[xs]['Water_lvl_elv'] = 100.01
            self.sim_water_lvl[xs]['Water_lvl_elv_high'] =  self.sim_water_lvl[xs]['Water_lvl_elv'] + 5.0
            self.sim_water_lvl[xs]['Water_lvl_elv_low'] = self.sim_water_lvl[xs]['Water_lvl_elv'] - 2.0
        
        
    
    def write_Testcd2(self):
        
        Testcd2 = open(self.stochICE.prjDir + '/' + 'DOUT7' + '/' + 'TESTCD2' + '.txt','w')
        
        def string_length_adjustment(string,desired_length):
            
            while len(string) < desired_length:
                
                string = ' ' + string
                
            return string
                                
                
        
        def write_Testcd2_reach_header(Reach_number):
            
            water_lvl_max = -9999.0
            water_lvl_min = 9999.0
            Manning = 0.0
            Reach_length = -9999.0
            Interp_interval = str(self.interInt) + '.'
            
            for xs_number in self.riv_xs_data:
                
                if self.riv_xs_data[xs_number]['Reach'] == Reach_number:
                    
                    if self.sim_water_lvl[self.riv_xs_data[xs_number]['Hecras xs']]['Water_lvl_elv_high'] > water_lvl_max:

                        water_lvl_max = self.sim_water_lvl[self.riv_xs_data[xs_number]['Hecras xs']]['Water_lvl_elv_high'] 
                        
                    if self.sim_water_lvl[self.riv_xs_data[xs_number]['Hecras xs']]['Water_lvl_elv_low'] < water_lvl_min:
                        
                        water_lvl_min = self.sim_water_lvl[self.riv_xs_data[xs_number]['Hecras xs']]['Water_lvl_elv_low']
                        
                    if int(self.riv_xs_data[xs_number]['Rivice Chainage']) > Reach_length : 
                        
                        Reach_length = self.riv_xs_data[xs_number]['Rivice Chainage']
                        
                        Manning = self.riv_xs_data[xs_number]['Manning']
                  
                             
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
            
            numbering = 1
            
            for xs_number in self.riv_xs_data:
                
                if self.riv_xs_data[xs_number]['Reach'] == Reach_number:
                                        
                    Reach_xs_number = str(numbering) 
                    water_lvl_elv = str(self.sim_water_lvl[self.riv_xs_data[xs_number]['Hecras xs']]['Water_lvl_elv'])   
                    discharge = str(self.sim_water_lvl[self.riv_xs_data[xs_number]['Hecras xs']]['Discharge'])
                    water_lvl_elv_high = str(self.sim_water_lvl[self.riv_xs_data[xs_number]['Hecras xs']]['Water_lvl_elv_high'])
                    water_lvl_elv_low = str(self.sim_water_lvl[self.riv_xs_data[xs_number]['Hecras xs']]['Water_lvl_elv_low'])
 
                    Reach_xs_number = string_length_adjustment(Reach_xs_number,9)
                    water_lvl_elv = string_length_adjustment(water_lvl_elv,11)
                    discharge = string_length_adjustment(discharge,10)
                    water_lvl_elv_high = string_length_adjustment(water_lvl_elv_high,8)
                    water_lvl_elv_low = string_length_adjustment(water_lvl_elv_low,8)
                    
                    Testcd2.write(Reach_xs_number + water_lvl_elv + discharge + '                                ' + water_lvl_elv_low + water_lvl_elv_high)
                    Testcd2.write('\n')
             
                numbering += 1
            
            
            if Reach_number < Number_of_reaches:
                
                Testcd2.write('\n')
                
            else:
                
                Testcd2.write('0')
                
        
        Number_of_reaches = 0
        
        for xs_number in self.riv_xs_data:
            
            if self.riv_xs_data[xs_number]['Reach'] > Number_of_reaches:
                
                Number_of_reaches = self.riv_xs_data[xs_number]['Reach']
    
    
        for Reach in range(1,Number_of_reaches + 1):
            
            write_Testcd2_reach_header(Reach)
            write_Testcd2_reach_data(Reach,Number_of_reaches)
            
                           
        Testcd2.close()                
                    
                
                
    def launch_Cd1xebat_file(self):
        
        print('')
        
        os.startfile(self.stochICE.prjDir + '/' + 'DOUT7' + '/' + 'Cd1x_e.bat')
        
        # os.startfile(self.stochICE.prjDir + '/' + 'DOUT7' + '/' + 'Cd1x_e.bat')
        
    def launch_Cd2pgmaexe_file(self):
        print('')
            
            
        
    def write_TAPE5(self):
        print('')
        
            
            
""" 
Notes

To create Dout7 -> Need Testcd2 &


"""        
            
            
            
            
            
            
            
            
            
            
            
            