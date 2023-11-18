"""
This is a testing file for stochICE v0.1
-It is easiest to place this file in the HEC-RAS project folder
"""


import os
import stochICE as ice


"""
User input parameters (eventually these will be selectable in a PyQt5 graphical user interface)
"""

path = os.getcwd() + "\Secteur_NeufPasMC_testing"
batch_ID="Test_2"
ras = "Secteur_neufpas.prj"
geo = "Secteur_neufpas.g01"
flowFile = "Secteur_neufpas.f03"
wse= 'WSE (13 avril 2011).Terrain.MNT_Point_a_neuf_pas.tif'

NSims = 10

thick=[0.3,0.7]
phi=[40,50]
flows=[100,300]

# locations=[[3731,520],[8370,5026],[3327,520]]
locations=[[8370,520]]




Roger=ice.stochICE(prjDir=path,
                                  batch_ID=batch_ID,
                                  ras_file=ras,
                                  geo_file=geo,
                                  flow_file=flowFile,
                                  wse_file=wse,
                                  NSims=NSims,
                                  thick_range=thick,
                                  phi_range=phi,
                                  flow_range=flows,
                                  locations=locations,
                                  code='RIVICE',
                                  clrRes=True,
                                  compRes=True,
                                  fun_mode=True)









    
# riv_xs_data = {}

# xs_number = 1
# xs_prec = ''

# if Roger.bridge:                        
    
#     bridge_number = 1
    
#     for xs in Roger.xs_data:
        
#         if xs_prec == '':
        
#             riv_xs_data[str(xs_number)] = {}
#             riv_xs_data[str(xs_number)]['Hecras xs'] = xs
            
            
#         elif bridge_number <= len(Roger.bridge_data) and Roger.bridge_data[str(bridge_number)]['chainage'] > float(xs):
            
#             riv_xs_data[str(xs_number)] = {}
#             riv_xs_data[str(xs_number)]['Hecras xs'] = xs_prec
            
#             xs_number += 1
            
#             riv_xs_data[str(xs_number)] = {}
#             riv_xs_data[str(xs_number)]['Hecras xs'] = xs
            
#             xs_number += 1
            
#             riv_xs_data[str(xs_number)] = {}
#             riv_xs_data[str(xs_number)]['Hecras xs'] = xs
            
#             bridge_number += 1
                   
            
#         elif Roger.xs_data[xs]['Manning']['val_MAIN'] != Roger.xs_data[xs_prec]['Manning']['val_MAIN'] :
            
#             riv_xs_data[str(xs_number)] = {}
#             riv_xs_data[str(xs_number)]['Hecras xs'] = xs_prec
            
#             xs_number += 1
            
#             riv_xs_data[str(xs_number)] = {}
#             riv_xs_data[str(xs_number)]['Hecras xs'] = xs
        
#         else :
            
#             riv_xs_data[str(xs_number)] = {}
#             riv_xs_data[str(xs_number)]['Hecras xs'] = xs
    
#         xs_number += 1        
#         xs_prec = xs
        
# else:
    
#     for xs in Roger.xs_data:
        
#         if xs_prec == '':
        
#             riv_xs_data[str(xs_number)] = {}
#             riv_xs_data[str(xs_number)]['Hecras xs'] = xs
            
                   
#         elif Roger.xs_data[xs]['Manning']['val_MAIN'] != Roger.xs_data[xs_prec]['Manning']['val_MAIN'] :
            
#             riv_xs_data[str(xs_number)] = {}
#             riv_xs_data[str(xs_number)]['Hecras xs'] = xs_prec
            
#             xs_number += 1
            
#             riv_xs_data[str(xs_number)] = {}
#             riv_xs_data[str(xs_number)]['Hecras xs'] = xs
        
#         else :
            
#             riv_xs_data[str(xs_number)] = {}
#             riv_xs_data[str(xs_number)]['Hecras xs'] = xs
    
#         xs_number += 1        
#         xs_prec = xs













# Roger.stochHECRAS.make_ensemble_flood_map()






# xsData=Roger.xsData
# BridgeData = Roger.BridgeData
# RivicexsData = Roger.stochRIVICE.riv_xsData


# bob.simInputPars
# bob.simProfileData


# TopIceMaxDepth=[]

# for key, value in bob.simProfileData.items():
#     print(TopIceMaxDepth.append(value['TopIceMaxDepth']))


# import numpy as np
# TopIceMaxDepth=np.asarray(TopIceMaxDepth)

# TopIceMaxDepth.mean(axis=0)
# TopIceMaxDepth.std(axis=0)
# TopIceMaxDepth.min(axis=0)
# TopIceMaxDepth.argmax(axis=0)

# TopIceMaxDepth.max(axis=0)

# TopIceMaxDepth[:,0]

# TopIceMaxDepth.mean(axis=0)[0]

# xs=10
# diff=TopIceMaxDepth[:,xs] - TopIceMaxDepth.mean(axis=0)[xs]


# # plt.scatter(diff,init_thick)


# discharge=np.asarray(bob.simInputPars['Q'])
# phi=np.asarray(bob.simInputPars['phi'])
# init_thick=np.asarray(bob.simInputPars['init_thick'])

# import pandas as pd
# import matplotlib.pyplot as plt
# data=pd.DataFrame()
# data["Q"]=discharge
# data["phi"]=phi
# data["init_thick"]=init_thick
# data["diff"]=diff


# corr = data.corr()
# fig = plt.figure()
# ax = fig.add_subplot(111)
# cax = ax.matshow(corr,cmap='coolwarm', vmin=-1, vmax=1)
# fig.colorbar(cax)
# ticks = np.arange(0,len(data.columns),1)
# ax.set_xticks(ticks)
# plt.xticks(rotation=90)
# ax.set_yticks(ticks)
# ax.set_xticklabels(data.columns)
# ax.set_yticklabels(data.columns)
# plt.show()




# np.corrcoef(diff,discharge)
# np.corrcoef(diff,phi)
# np.corrcoef(diff,init_thick)



# import matplotlib.pyplot as plt



# _ = plt.hist(hist)  # arguments are passed to np.histogram

# plt.title("Histogram with 'auto' bins")
# Text(0.5, 1.0, "Histogram with 'auto' bins")

# plt.show()


# bob.simProfileData['sim_1']['WSE']-bob.simProfileData['sim_1']['MinChEle']

# bob.simProfileData['sim_2']['WSE']-bob.simProfileData['sim_2']['MinChEle']


# for item in xsData.items():
#     print(xsData[item]["chainage"])
# dataIneed=xsData['8054']

# class testClass2000():

#     def __init__(self):
        
#         print('This is a simple class')
#         bobsName="Bob"

#     def thisFunction(self):
        
#         print("This function")
        
        
        

# bob=testClass2000()
# bob.bobsName

# bob.bobsname="Paul"

# bob.bobsname
# bob.thisFunction()
# simData=bob.simInputPars
# simProfileData=bob.simProfileData


# bob.chainages
# import matplotlib.pyplot as plt
# import numpy as np

# plt.plot(bob.chainages,simProfileData['sim_4']['WSE'])

# mean_wse=[]
# min_wse=[]
# max_wse=[]
# for key in simProfileData.keys():
#     print(key)

#     mean_wse.append(np.mean(simProfileData[key]['WSE']))
#     min_wse.append(min(simProfileData[key]['WSE']))
#     max_wse.append(max(simProfileData[key]['WSE']))

# import numpy as np

# StrRS=bob.RC.Geometry_GetNodes(bob.RiverID,bob.ReachID)[3]

# bob.createStocasticFloodRaster()

# dir(bob.RC)

# Section = np.empty([bob.NNod],dtype=float)
# for n in range(0,bob.NNod):
#     Section[n] = float(StrRS[n])






# [k for k,v in simData.items() print(k) >= 17]

# wse=bob.TabWSE






"""
Extract HEC-RAS data at nodes
"""

# StrRS = RC.Geometry_GetNodes(RiverID, ReachID)[3];
# Section = np.empty([NNod],dtype=float)
# for n in range(0,NNod):
#   Section[n] = float(StrRS[n])


# """
# Save simulation parameters to file
# """

# par_path=project_directory+"\SimulationParameters"

# if not os.path.exists(par_path):
#     os.makedirs(par_path)

# arr_input = np.asarray([TabFLOW, TabICEInput,TabFriction])
# np.savetxt(par_path+"\inputMC.csv", arr_input, delimiter=",",fmt="%3.2f")
# np.savetxt(par_path+"\TabWSE.csv", TabWSE, delimiter=",",fmt="%3.2f")
# np.savetxt(par_path+"\TabICE.csv", TabWSE, delimiter=",",fmt="%3.2f")
# np.savetxt(par_path+"\Section.csv", Section, delimiter=",",fmt="%d")


# """
# -----------------------------------------------------------
# -----------------------------------------------------------
# Postprocess results
# -----------------------------------------------------------
# -----------------------------------------------------------
# """


# """
#               ``'-.,_,.-'``'-.,_,.='``'-.,_,.-'``'-.,_,.='``
#               ``'-.,_,.-'``'-.,_,.='``'-.,_,.-'``'-.,_,.='``
#               ``'-.,_,.-'``'-.,_,.='``'-.,_,.-'``'-.,_,.='``
#               ----------Plot water level profiles-----------
#               ``'-.,_,.-'``'-.,_,.='``'-.,_,.-'``'-.,_,.='``
#               ``'-.,_,.-'``'-.,_,.='``'-.,_,.-'``'-.,_,.='``
#               ``'-.,_,.-'``'-.,_,.='``'-.,_,.-'``'-.,_,.='``
# """

#


# """
#                    Produce flood map
#                       _,--',   _._.--._____
#                .--.--';_'-.', ";_      _.,-'
#               .'--'.  _.'    {`'-;_ .-.>.'
#                     '-:_      )  / `' '=.
#                       ) >     {_/,     /~)
#                       |/               `^ .'
# """

# graphique stochastique




# """
# Statistical analysis of results
# """


# #Calculs des PDA

# tiff_stoch = rasterio.open(floodmap_output)
# arr_stoch = tiff_stoch.read()
# x = np.nonzero(arr_stoch > 1)

# all_im = np.zeros((1000,np.shape(x)[1]))
# listdir = os.listdir(tif_path)

# for count, item in enumerate(listdir):
#       data_name = tif_path + "\\"+ item
#       tiff = rasterio.open(data_name)
#       arr2 = tiff.read()
#

# for m in range(4,len(listdir)):
#       data_name = tif_path + "\\"+ listdir[m]
#       tiff = rasterio.open(data_name)
#       arr2 = tiff.read()
#
#       for mm in range(1,np.shape(x)[1]):
#               all_im[m,mm]=arr2[0,x[1][mm],x[2][mm]]

# all_im[all_im<1]=np.nan

# WSE_100 = np.zeros((np.shape(x)[1]))
# WSE_20 = np.zeros((np.shape(x)[1]))
# WSE_2 = np.zeros((np.shape(x)[1]))

# for mm in range(1,np.shape(x)[1]):
#   WSE_100[mm] = np.nanquantile(all_im[:,mm], 0.99)
#   WSE_20[mm] = np.nanquantile(all_im[:,mm], 0.95)
#   WSE_2[mm] = np.nanquantile(all_im[:,mm], 0.50)


# arr_100 = np.zeros((1,np.shape(arr2)[1],np.shape(arr2)[2]))
# arr_20 = np.zeros((1,np.shape(arr2)[1],np.shape(arr2)[2]))
# arr_2 = np.zeros((1,np.shape(arr2)[1],np.shape(arr2)[2]))


# for mm in range(1,np.shape(x)[1]):
#   arr_100[0,x[1][mm],x[2][mm]]=WSE_100[mm]
#   arr_20[0,x[1][mm],x[2][mm]]=WSE_20[mm]
#   arr_2[0,x[1][mm],x[2][mm]]=WSE_2[mm]

# data_out = r"C:\Users\dugj2403\Desktop\Secteur_NeufPasMC\MonteCarlo/WSE_100.tif"
# result = rasterio.open(data_out,'w',driver='GTiff',height=tiff.shape[0],width=tiff.shape[1],count=1,dtype=arr_100.dtype,crs=tiff.crs,transform=tiff.transform)
# result.write(arr_100[0,:,:],1)
# result.close()

# data_out = r"C:\Users\dugj2403\Desktop\Secteur_NeufPasMC\MonteCarlo/WSE_20.tif"
# result = rasterio.open(data_out,'w',driver='GTiff',height=tiff.shape[0],width=tiff.shape[1],count=1,dtype=arr_100.dtype,crs=tiff.crs,transform=tiff.transform)
# result.write(arr_20[0,:,:],1)
# result.close()

# data_out = r"C:\Users\dugj2403\Desktop\Secteur_NeufPasMC\MonteCarlo/WSE_2.tif"
# result = rasterio.open(data_out,'w',driver='GTiff',height=tiff.shape[0],width=tiff.shape[1],count=1,dtype=arr_100.dtype,crs=tiff.crs,transform=tiff.transform)
# result.write(arr_2[0,:,:],1)
# result.close()
