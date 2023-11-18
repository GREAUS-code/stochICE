# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 13:52:46 2023

@author: dugj2403
"""

"TEST DE GITHUB"

class StochRIVICE():
    
    # pass-in stochICE instance 
    def __init__(self,stochICE,interInt):
        
        self.stochICE=stochICE
        self.interInt=interInt
           
    def get_riv_xsections(self): 
        
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
            
            
    def RiviceChainage(self):
        
        def mon_arrondi(x, base):
            return base * round(x/base)
        
        flag = False
        reste_arrondi = 0            
        xs_number_precedent = ''
        xs_number_precedent_precedent = ''
        
        for xs_number in self.riv_xs_data:
            
            if xs_number == '1':
                self.riv_xs_data[xs_number]['Rivice Chainage'] = 0
                
            elif xs_number == str(len(self.riv_xs_data)):
                self.riv_xs_data[xs_number]['Rivice Chainage'] = mon_arrondi(self.riv_xs_data[xs_number_precedent]['Rivice Chainage'] + (int(self.riv_xs_data[xs_number_precedent]['Hecras xs']) - int(self.riv_xs_data[xs_number]['Hecras xs'])), self.interInt)
                 
            elif flag:
                self.riv_xs_data[xs_number]['Rivice Chainage'] = self.riv_xs_data[xs_number_precedent]['Rivice Chainage'] + (int(self.riv_xs_data[xs_number_precedent]['Hecras xs']) - (reste_arrondi + int(self.riv_xs_data[xs_number]['Hecras xs'])))
                flag = False
                
            elif self.riv_xs_data[xs_number]['Hecras xs'] == self.riv_xs_data[xs_number_precedent]['Hecras xs'] :
                
                self.riv_xs_data[xs_number]['Rivice Chainage'] = mon_arrondi(self.riv_xs_data[xs_number_precedent_precedent]['Rivice Chainage'] + (int(self.riv_xs_data[xs_number_precedent_precedent]['Hecras xs']) - int(self.riv_xs_data[xs_number_precedent]['Hecras xs'])), self.interInt)
                
                self.riv_xs_data[xs_number_precedent]['Rivice Chainage'] = self.riv_xs_data[xs_number]['Rivice Chainage']   
                
                reste_arrondi = mon_arrondi(self.riv_xs_data[xs_number_precedent_precedent]['Rivice Chainage'] + (int(self.riv_xs_data[xs_number_precedent_precedent]['Hecras xs']) - int(self.riv_xs_data[xs_number_precedent]['Hecras xs'])), self.interInt) - (self.riv_xs_data[xs_number_precedent_precedent]['Rivice Chainage'] + (int(self.riv_xs_data[xs_number_precedent_precedent]['Hecras xs']) - int(self.riv_xs_data[xs_number_precedent]['Hecras xs'])))
                flag = True                                                                                    
            
            else:
                self.riv_xs_data[xs_number]['Rivice Chainage'] = self.riv_xs_data[xs_number_precedent]['Rivice Chainage'] + (int(self.riv_xs_data[xs_number_precedent]['Hecras xs']) - int(self.riv_xs_data[xs_number]['Hecras xs']))        
            
            xs_number_precedent_precedent = xs_number_precedent
            
            xs_number_precedent = xs_number
    
     
    def RiviceManning(self):
        
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
            
            
    def RiviceXSectionMainChannelGeometry(self):
        
        for xs_number in self.riv_xs_data:
            self.riv_xs_data[xs_number]['MainChannelGeometry'] = self.stochICE.xs_data[self.riv_xs_data[xs_number]['Hecras xs']]['MainChannelGeometry']['xy']
            
            
    def RiviceReach(self):
        
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
            
            
    def RiviceXSectionID(self):
        
        for xs_number in self.riv_xs_data:
            
            self.riv_xs_data[xs_number]['XS ID'] = int(xs_number)*1000
        
        
        
        
    def RiviceDistanceXSectionPrecedente(self):
        
        xs_number_precedent = ''
        
        for xs_number in self.riv_xs_data:
            
            if xs_number == '1':
            
                self.riv_xs_data[xs_number]['Distance XS Precedente'] = 0
            
            else:
                
                self.riv_xs_data[xs_number]['Distance XS Precedente'] = self.riv_xs_data[xs_number]['Rivice Chainage'] - self.riv_xs_data[xs_number_precedent]['Rivice Chainage']
                
                
            xs_number_precedent =xs_number   
                
                
   

    def write_Cd1test(self):    
        
        
        Cd1test = open(self.stochICE.prjDir + '/' + 'CD1TEST' + '.txt','w')
        
        
        def Cd1testReachHeader(xs_number):
                       
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
            
            
        def Cd1testXSectionHeader(xs_number):
            
            Header_XS_ID = str(self.riv_xs_data[xs_number]['XS ID']) + '.'
            
            while len(Header_XS_ID) < 10:
                
                Header_XS_ID = ' ' + Header_XS_ID
                
                
            Header_Distance_xs_prec = str(self.riv_xs_data[xs_number]['Distance XS Precedente']) + '.00'
                
            while len(Header_Distance_xs_prec) < 11 : 
                
                Header_Distance_xs_prec = ' ' + Header_Distance_xs_prec
                
            
            Header_XS = 'SECT' + Header_XS_ID + Header_Distance_xs_prec
              
            Cd1test.write(Header_XS)
            Cd1test.write('\n')
            
            
        def Cd1testXSectionGeometry(xs_number):
            
            
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
                
                Cd1testReachHeader(xs_number)
                Cd1testXSectionHeader(xs_number)
                Cd1testXSectionGeometry(xs_number)
                
            
            elif self.riv_xs_data[xs_number]['Reach'] != self.riv_xs_data[xs_number_precedent]['Reach']: 
                
                Cd1testReachHeader(xs_number)
                Cd1testXSectionHeader(xs_number)
                Cd1testXSectionGeometry(xs_number)           
            
            
            else:
                
                Cd1testXSectionHeader(xs_number)
                Cd1testXSectionGeometry(xs_number)
                
                        
            xs_number_precedent = xs_number