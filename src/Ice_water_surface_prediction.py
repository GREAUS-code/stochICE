from gekko import GEKKO
import math
import pandas as pd
import numpy as np

import stochICE as ice
import matplotlib.pyplot as plt

"""
Why does this not work?
"""

"""
Find downstream water level with ice cover for corresponding discharge 
selected at upstream boundary in RIVICE. 
"""

# Manning's solver Credited to https://github.com/alexiusacademia/ChannelFlowLib




"""
Get bathymetry for downstream cross-section. This will be taken from xs_data in
stochICE in the final version. 
"""

bathy='downstream_cross_section.txt'
# bathy='upstream_cross_section.txt'

def readBathy(path):
    
    df=pd.read_csv(path,delim_whitespace=True,names=['station','elevation'])
     
    return df

bathy_pnts=readBathy(bathy)
pts = tuple(bathy_pnts.loc[:, ['station', 'elevation']].to_records(index=False))


"""
Setup channel
"""

channel = ice.irregularSection.IrregularSection(pts)
channel.set_average_rougness(0.03)
channel.set_bed_slope(0.00045)
channel.set_water_elevation(67.7)
channel.analyze()
channel.discharge
channel.wetted_area
channel.wetted_perimeter

"""
Create a rating curve
"""
max_discharge=500
max_elev = bathy_pnts['elevation'].max()
min_elev = bathy_pnts['elevation'].min()
interval = 0.03
intervals = ((max_elev - min_elev) / interval)

elevs = []
discharges = []
wetted_perimeters=[]
wetted_areas=[]
top_widths=[]


for i in range(int(intervals) + 1):
    
    #seems like keeping 2 decimal places causes errors, so I only keep 1
    elev = round(min_elev + ((i+1) * interval),1)
    
    channel.set_water_elevation(elev)
    
    channel.analyze()
    print(channel.state)
    if bool(channel.state) == True:
        
        if channel.discharge <= max_discharge:
            elevs.append(elev)
            discharges.append(channel.discharge)
            wetted_perimeters.append(channel.wetted_perimeter)
            wetted_areas.append(channel.wetted_area)
            top_widths.append(channel.top_width)
        else:
            break
    
    if bool(channel.state) == False:
        break
        

plt.plot(discharges, elevs)
plt.show()


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



