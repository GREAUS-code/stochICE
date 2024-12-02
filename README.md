# stochICE

Python package allowing either HECRAS to stochastically model the effect of river ice jams on water levels using a Monte-Carlo approach. The floodmaps produced provide a sense of the probability that an area of interest will be flooded. The user conveniently setups the geometry of the hydraulic model using the RASMapper module available within HECRAS. 

## HECRAS - ice modelling overview

The wide-river ice jam functionality of HEC-RAS allows the user to predict the thickness of an ice jam that forms between user specified cross-sections and the effects the ice jam has on the river's water levels. The model calculates the thickness ($t$, m) of the jam at each section by balancing the forces preventing the ice jam from moving downstream with the external forces attempting to dislodge it. During model development in HECRAS, the user must specify the following: discharge, downstream water surface level (usually based on the normal depth in the reach), locations of the upstream and downstream end of the ice jam, open water sections, sections with an ice cover, manning coefficients of bed and ice, ice cover thickness, ice jam porosity (as a void fraction), the ice jam's internal friction angle and the maximum flow velocity under the jam.

## Required inputs:

- the most recent version of HECRAS installed on a Windows operating system
- a digital elevation model of the study reach with bathymetry embedded (e.g. provincial LiDAR dataset)
- discharge and ice cover thickness distributions and locations of potential ice jam location (among other possible parameter distributions)

## Outputs include: 

- longitudinal water surface profiles
- longitudinal ice thickness profiles  
- 2D water-level distribution geotiffs of individual ice jam events
- 2D statistical water-level and maximum depths over the ensemble of simulated ice jams.
![(see example here)](https://github.com/GREAUS-code/stochICE/tree/main/imgs/stochMapExample.jpg?raw=true)

- descriptive statistics of the ensemble of tested ice jam scenarios

