"""
This is a tutorial/testing input script for running RIVICE stochastically with stochICE
"""

import sys
import os
import stochICE as ice
import pandas as pd
import math



"""
Step 1: setup the environment
"""
# stochICE_path = r"C:\Users\jason\Desktop\stochICE\src"
# sys.path.append(os.path.abspath(stochICE_path)) #adds stochICE to path
# import stochICE as ice



"""
Step 2: specify input files and paths
"""


#Change to the the path containing your HECRAS files
path = r'C:\Users\jason\Desktop\stochICE\examples\tut_1_single_reach'
batch_ID="Secteur_neufpas"
ras = "Secteur_neufpas.prj"
geo = "Secteur_neufpas.g01"
flowFile = "Secteur_neufpas.f01"
wse= 'WSE (13 avril 2011).Terrain.MNT_Point_a_neuf_pas.tif'



"""
Step 3: specify RIVICE simulation parameters
"""

#Number of simulations 
NSims = 10

#Simulation batch size (NSims/GrpSize must not have a remainder)
GrpSize=10

#Inputs necessary for HECRAS precursor openwater simulation
thick=[0,0]
phi=[0,0]
flows=[70,70]            # best to choose a discharge that is on the larger side of the distribution 
locations=[[0,0]]        # 0's remove ice jams from HECRAS precursor simulation
ds_slope=0.00031
max_Q=350


#RIVICE simulation parameters
riv_sim_days = 2         # days
riv_timestep = 30        # seconds
riv_ice_start=0.5        # days 
riv_ice_end=2            # days
riv_interInt = 30        # x-section interpolation interval (m) (recommendation: width of the river)
riv_profile_interval=0.5 #Profile time output interval in days (must be a divisor of riv_sim_days)
riv_dwn_bc_opt = 1       # 1 - Normal depth, adjusted for ice cover, 2 - as a stochastic variable (define distribution in stoch_variables dictionary)



# #RIVICE simulation parameters
# riv_sim_days = 0.208         # days
# riv_timestep = 1        # seconds
# riv_ice_start=0.083        # days 
# riv_ice_end=0.208            # days
# riv_interInt = 20        # x-section interpolation interval (m) (recommendation: width of the river)
# riv_profile_interval=0.208 #Profile time output interval in days (must be a divisor of riv_sim_days)
# riv_dwn_bc_opt = 1       # 1 - Normal depth, adjusted for ice cover, 2 - as a stochastic variable (define distribution in stoch_variables dictionary)


"""
Step 4: specify variables to stochastically model. 
- Currently, variables can only be selected from a normal distribution. 
- Future work will include other distributions. 
"""

#Syntax: 'VariableName in TAPE5.txt':[mean, std, nmb_samples], see init_default_ice_parms() in stochRIVICE_JD for a complete list.
#At the moment, only a normal distribution has been implemented. 
# stoch_variables={'Frontthick':[0.3,0.1,5000],
#                  'Q':[50,10,5000],
#                  'IceVol':[270,30,5000],
#                  'RLOCBRG':[7750,3000,5000],
#                  'DAYSBR':[0.5,0,5000], 
#                  'FACTOR3':[0.025,0.01,5000]}

stoch_variables={'Frontthick':['uniform',[0.1,0.6],2],
                 'Q':['uniform',[50,100],2],
                 'IceVol':['uniform',[1,10],2],
                 'RLOCBRG':['uniform',[1000,8000],0],
                 'FACTOR3':['uniform',[0.01,0.08],2]
                 }

"""
Step 5: Setup up simulation folders and write TAPE5.txt for each simulation.
"""
Roger=ice.stochICE(prjDir=path,

                                  batch_ID=batch_ID,
                                  ras_file=ras,
                                  geo_file=geo,
                                  flow_file=flowFile,
                                  wse_file=wse,
                                  NSims=NSims,
                                  GrpSize=GrpSize,
                                  thick_range=thick,
                                  phi_range=phi,
                                  flow_range=flows,
                                  ds_slope=ds_slope,
                                  max_Q=max_Q,
                                  locations=locations,
                                  riv_dwn_bc_opt=riv_dwn_bc_opt,
                                  code='RIVICE',
                                  days=riv_sim_days,
                                  timestep=riv_timestep,
                                  ice_start=riv_ice_start,
                                  ice_end=riv_ice_end,
                                  interval=riv_interInt,
                                  profile_int=riv_profile_interval,
                                  clrRes=True,
                                  compRes=True,
                                  fun_mode=False,
                                  sleep=5,
                                  stochvars=stoch_variables)



"""
Step 6: Launch simulations in parallel
"""

Roger.stochRIVICE.launch_RIVICE_parallel()


"""
Step 7: Plot water surface profile envelope for reach
"""

Roger.stochRIVICE.plot_profiles()

profiles=Roger.stochRIVICE.sim_profiles

"""
Step 8: Plot exceedance probability for a specified chainage along reach
"""

chainage=4000
Roger.stochRIVICE.plot_prob_exceedance(chainage)




"""
Code for general distribution spec
"""
import random


                

def 
variable_values={}

for variable, parms in test_variables.items():
    
    if parms[0] == 'uniform':

        secure_random = random.SystemRandom()
        random_value = round(secure_random.uniform(parms[1][0], parms[1][1]), parms[2])
        variable_values[variable]=random_value
        
    if parms[0] == 'normal':
        
        """
        Careful, this code can produce unphysical negative values!
        """

        rng = np.random.default_rng()
        values = rng.normal(parms[1][0], parms[1][1], size=parms[1][2])
        variable_values[variable]=random.choice(values)

    if parms[0] == 'Gumbel': 
        
        print('Gumbel distributions not yet implemented')
       

plt.hist(values)




    
    

    



    
    




































