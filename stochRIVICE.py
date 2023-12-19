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
                
                
    def get_reach_data(self):
        
        self.reach_data = {}
        
        #Recueil du numero de reach et du chainage de sa limite aval
        
        
        for i in range(1,len(self.riv_xs_data) + 1):
            
            if i == len(self.riv_xs_data):
                
                Reach_number = self.riv_xs_data[str(i)]['Reach']
                Reach_end_chainage = self.riv_xs_data[str(i)]['Rivice Chainage']
                
                self.reach_data[str(Reach_number)] = {}
                self.reach_data[str(Reach_number)]['Reach_number'] = Reach_number
                self.reach_data[str(Reach_number)]['Reach_end_chainage'] = Reach_end_chainage
                
                
            else:
            
                Reach_number = self.riv_xs_data[str(i)]['Reach']
                Reach_end_chainage = self.riv_xs_data[str(i)]['Rivice Chainage']
                
                Next_reach_number = self.riv_xs_data[str(i+1)]['Reach']
                
                
                if Reach_number != Next_reach_number:
                    
                    self.reach_data[str(Reach_number)] = {}
                    self.reach_data[str(Reach_number)]['Reach_number'] = Reach_number
                    self.reach_data[str(Reach_number)]['Reach_end_chainage'] = Reach_end_chainage
                    
        #Recueil de la longueur de chaque reach et du nombre de XS qu'ils comprennent
        #incluant les XS interpolees            
        for i in range(1,len(self.reach_data)+1):
            
            if i == 1:
                
                Reach_length = self.reach_data[str(i)]['Reach_end_chainage']                
                Num_of_XS = int(Reach_length)/self.interInt + 1
                
                self.reach_data[str(i)]['Reach_length'] = Reach_length                
                self.reach_data[str(i)]['Number_of_XS'] = Num_of_XS
                
            else:
                            
                Reach_end_chainage = self.reach_data[str(i)]['Reach_end_chainage']
                Prev_reach_end_chainage = self.reach_data[str(i-1)]['Reach_end_chainage']               
                
                Reach_length = Reach_end_chainage - Prev_reach_end_chainage                
                Num_of_XS = Reach_length/self.interInt + 1
                
                self.reach_data[str(i)]['Reach_length'] = Reach_length                               
                self.reach_data[str(i)]['Number_of_XS'] = Num_of_XS
                
                
        

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
            
            def water_lvl_format_adjusment(string):
                
                if len(string) == 4:
                    
                    string = string + '0'
                    
                return string
            
            numbering = 1
            
            for xs_number in self.riv_xs_data:
                
                if self.riv_xs_data[xs_number]['Reach'] == Reach_number:
                                        
                    Reach_xs_number = str(numbering) 
                    water_lvl_elv = str(self.sim_water_lvl[self.riv_xs_data[xs_number]['Hecras xs']]['Water_lvl_elv'])   
                    discharge = str(self.sim_water_lvl[self.riv_xs_data[xs_number]['Hecras xs']]['Discharge'])
                    water_lvl_elv_high = str(self.sim_water_lvl[self.riv_xs_data[xs_number]['Hecras xs']]['Water_lvl_elv_high'])
                    water_lvl_elv_low = str(self.sim_water_lvl[self.riv_xs_data[xs_number]['Hecras xs']]['Water_lvl_elv_low'])
                    
                    water_lvl_elv = water_lvl_format_adjusment(water_lvl_elv)
                    water_lvl_elv_high = water_lvl_format_adjusment(water_lvl_elv_high)
                    water_lvl_elv_low = water_lvl_format_adjusment(water_lvl_elv_low)
 
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
                    
                
                
    def launch_Cd1xebat_Cd2pgmaexe_file(self):
        
        print('')
        
        os.startfile(self.stochICE.prjDir + '/' + 'DOUT7' + '/' + 'Cd1x_e_&_Cd2pgm_a.bat')            
            
        
    def write_TAPE5(self):
        
        def string_length_adjustment(string,length,position):
            
            if position == "F":
            
                while len(string) < length:
                    
                    string = " " + string
            
            elif position == "B":
            
                while len(string) < length:
                    
                    string = string + " "
                      
            else:               
                print('Error : The position (3rd argument), should be F (for front adjustment) or B (for back adjustment)')
                
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
        NINC = 2880
        PERIOD = 86400
        RATIO = ""
        MAXITR = "" 
        EPS = ""
        LPER = ""
        NRINC = ""
        
        # Network Parameters (voir le manuel de l'utilisateur de RIVICE pour plus de détails sur ces options)       
        NREACH = self.riv_xs_data[str(len(self.riv_xs_data))]["Reach"] # -> Cette variable n'a pas besoin d'être précisée par l'utilisateur
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
        Q_inputs = [[0.0,2880.0],\
                    [100.0,100.0]]            
        Q_Type = 2 #Boundary condition type
        Q_Time_dep = 2 #Boundary condition time dependance
        Q_Serie_len = len(Q_inputs[0]) #Boundary condition serie length
        Q_Intr_type = 1 #Interpolation type used to interpolate between the
                        #boundary condition serie terms              
        
        
        # Downstream water surface elevation boundary condition parameters
        H_inputs = [[0.0,2880.0],\
                    [66.08,66.08]]
        H_Type = 1 #Boundary condition type
        H_Time_dep = 2 #Boundary condition time dependance
        H_Serie_len = len(Q_inputs[0]) #Boundary condition serie length
        H_Intr_type = 1 #Interpolation type used to interpolate between the
                        #boundary condition serie terms   
        
        # Hydraulic hydrographs and profiles parameters
        HG_xs_chainage = [0.0,1685.0,4090.0,7095.0,8285.0]
        PF_time_step = [0.0,1440.0,2880.0] # Corresponds to the time step number
                                           # where the profiles will be
                                           # printed, not the real time
                                           

        
        
        """
        Écriture du fichier de contrôle TAPE5
        
        """
        
        TAPE5 = open(self.stochICE.prjDir + '/TAPE5.txt','w')
         
        # Writing project name
        if len(Project_title) > 80:
            
            print("Error : The project title should be 80 characters or less. Please shorten the project title")
        
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
        PERIOD = string_length_adjustment(str(PERIOD) + ".",10,"F")
            
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
        Dout7 = open(self.stochICE.prjDir + '/' + 'DOUT7' + '/' + 'DOUT7.txt','r')
        
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
                Reach_xs_num = string_length_adjustment("-" + str(Reach_xs_num),10,'F')
                
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
        US_Q_boundary_cdn = [Q_inputs,Q_Type,Q_Time_dep,Q_Serie_len,Q_Intr_type] #US boundary condition parameters
        
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
            
        i += 1
            
        reach_num = string_length_adjustment(str(i),10,'F')
        
        line = reach_num
        
        # Downstream boundary condition definition
        DS_WSE_boundary_cdn = [H_inputs,H_Type,H_Time_dep,H_Serie_len,H_Intr_type]
        
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
            
            for j in range(1,len(self.riv_xs_data)+1):
                
                if int(HG_xs_chainage[i]) == self.riv_xs_data[str(j)]['Rivice Chainage']:
                    
                    reach_num = self.riv_xs_data[str(j)]['Reach']                
                    xs_chainage = HG_xs_chainage[i]
                    
                    reach_num = string_length_adjustment(str(reach_num),10,'F')
                    xs_chainage = string_length_adjustment(str(xs_chainage),10,'F')
                    line = "         1         1         1"                    
                    
                    TAPE5.write(reach_num + xs_chainage + line)
                    TAPE5.write("\n")
                    
                    break
                    
        # Profiles parameters
        PF_num = len(PF_time_step)
        Reach_num = len(self.reach_data)
        
        TAPE5.write(string_length_adjustment(str(PF_num*Reach_num),10,'F'))
        TAPE5.write("\n")
        
        for i in range(PF_num):
            
            time_step = string_length_adjustment(str(PF_time_step[i]),10,'F')
            
            for j in range(1,Reach_num+1):
                
                Reach = string_length_adjustment(str(j),10,'F')               
                line = "         1         1         1"
                
                TAPE5.write(Reach + time_step + line)
                TAPE5.write("\n")
                
                
                
        # Writing Water quality boundary conditions
        TAPE5.write("WATER QUALITY BOUNDARY CONDITIONS")
        TAPE5.write("\n")
        
        
        # Remettre les lignes de codes commentées dans le blocs de variables
        # qui doivent être spécifiée par l'utilisateur et continuer le code
        # # Upstream water temperature boundary condition parameters
        # # 0 - Time step number
        # # 1 - Water temperature (Fahrenheit)
        # # 2 - Dispersive flux
        # # 3 - Total flux       
        # T_inputs = [[0.0,1440.0,2880.0],\
        #             [32.0,32.0,32.0],\
        #             [0.0,0.0,0.0],\
        #             [0.0,0.0,0.0]]
        # T_Type = 1 #Boundary condition type
        # T_Time_dep = 2 #Boundary condition time dependance
        # T_Serie_len = len(T_inputs[0]) #Boundary condition serie length
        
        
        
                
                
        # Writing Water quality time graphs and profiles
        TAPE5.write("WATER QUALITY TIME GRAPHS & PROFILES")
        TAPE5.write("\n")
        
        # Time graphs parameters
        
        
        
        
        # Profiles parameters
                
                
                
                
            
                
                    
                    
        
          
            

            
        
        
        
        
        
        
        


        
        
            
           
        
        
        
        
        
        
        
        
        
        

        
        TAPE5.close()
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
            
            
            
            
            
            
            
            
            
            
            
            