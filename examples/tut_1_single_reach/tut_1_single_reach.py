"""
This is a tutorial/testing input script for running RIVICE stochastically with stochICE
"""

import sys
import os
import stochICE as ice
import pandas as pd


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
flowFile = "Secteur_neufpas.f03"
wse= 'WSE (13 avril 2011).Terrain.MNT_Point_a_neuf_pas.tif'



"""
Step 3: specify RIVICE simulation parameters
"""

#Number of simulations 
NSims = 1

#Simulation batch size (NSims/GrpSize must not have a remainder)
GrpSize=1

#Inputs necessary for HECRAS precursor openwater simulation
thick=[0,0]
phi=[0,0]
flows=[80,80]            # best to choose a discharge that is on the larger side of the distribution 
locations=[[0,0]]        # 0's remove ice jams from HECRAS precursor simulation
ds_slope=0.00031
max_Q=350

#RIVICE simulation parameters
riv_sim_days = 2         # days
riv_timestep = 30        # seconds
riv_ice_start=0.5        # days 
riv_ice_end=2            # days
riv_interInt = 50        # x-section interpolation interval (m) (recommendation: width of the river)
riv_profile_interval=0.5 #Profile time output interval in days (must be a divisor of riv_sim_days)



"""
Step 4: specify variables to stochastically model. 
- Currently, variables can only be selected from a normal distribution. 
- Future work will include other distributions. 
"""

#Syntax: 'VariableName in TAPE5.txt':[mean, std, nmb_samples], see init_default_ice_parms() in stochRIVICE_JD for a complete list.
#At the moment, only a normal distribution has been implemented. 
stoch_variables={'Frontthick':[0.2,0.0,5000],'Q':[80,30,5000],'DWSE':[65.4,0,5000],'IceVol':[270,30,5000],'RLOCBRG':[130,0,5000],'DAYSBR':[0.5,0,5000]}


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


Roger.stochRIVICE.get_downstream_xs_data()
Roger.stochRIVICE.make_downstream_stage_discharge_curve()



import scipy.interpolate


"""
Add ice cover
"""

wp_ice=list(np.asarray(wetted_perimeters)+np.asarray(top_widths))


"""
Input variables
"""
n=0.03
n_ice=0.04
s=0.000450 # need this from hecras model (or ask the user to manually input it for each case)
t=0.2
target_Q=150

n_c=((1+(n_ice/n)**(3/2))/2)**(2/3)*n


ice_cover_discharges=[]
for i, area in enumerate(wetted_areas):
    
    V=(1/n_c)*((area/wp_ice[i])**(2/3))*(s)**0.5
    Q=V*area
    print(Q)
    ice_cover_discharges.append(Q)

y_interp = scipy.interpolate.interp1d(ice_cover_discharges, elevs)

print(y_interp(target_Q)+t*0.91)


y_interp(target_Q)
"""
Rectangular cross-section approximation
"""

# Define initial guesses
D_i=0.01
A_glace=B_rivice*D_i
P_glace=2*B_rivice+2*D_i
R_glace=A_glace/P_glace
V_glace=1/(n_c)*(R_glace**(2/3))*(s**(1/2))

# setup up equation solver
m = GEKKO()  

D_i = m.Var(value = D_i, lb = 0)
A_glace = m.Var(value = A_glace, lb = 0)
P_glace = m.Var(value = P_glace)
R_glace = m.Var(value = R_glace)
V_glace = m.Var(value = V_glace)

# solve equations
m.Equation(B_rivice*D_i==A_glace)
m.Equation(2*B_rivice+2*D_i==P_glace)
m.Equation(A_glace/P_glace == R_glace)
m.Equation(1/(n_c)*(R_glace**(2/3))*(s**(1/2)) == V_glace)
m.Equation(V_glace*A_glace == Q_rivice)

m.solve(disp = True)

wse_ice=WSE_rivice-(A_rivice/B_rivice)+D_i.value[0]+t
print(D_i.value[0], V_glace.value[0]*A_glace.value[0], wse_ice)





"""
Step 6: Launch simulations in parallel
"""

Roger.stochRIVICE.launch_RIVICE_parallel()


"""
Step 7: Plot water surface profile envelope for reach
"""

Roger.stochRIVICE.plot_profiles()

"""
Step 8: Plot exceedance probability for a specified chainage along reach
"""

chainage=3000
Roger.stochRIVICE.plot_prob_exceedance(chainage)

Roger.













    
    

    



    
    




































