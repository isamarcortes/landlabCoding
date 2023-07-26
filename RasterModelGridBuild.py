#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 11:25:46 2023

@author: isamarcortes
"""

from landlab import RasterModelGrid, imshow_grid
from landlab.components import LinearDiffuser
import numpy as np
import matplotlib.pyplot as plt
import rasterio as rio





File = rio.open('/Users/isamarcortes/Dropbox/Isamar/Papers_In_Prep/Paper_4/RasterLayersForLandlab/PR1.tif')
veg_array = File.read(1)
test = File.read(1)


sal_array = np.select([veg_array == 0, veg_array == 1], [35, 36], veg_array)



Island = RasterModelGrid((veg_array.shape))
vegetation = Island.add_field('vegetation',veg_array,at='node')



#salinity = Island.add_field('salinity',sal_array,at='node')
#Island.status_at_node[salinity == 35]=Island.BC_NODE_IS_FIXED_VALUE
#Island.status_at_node[salinity==36]=Island.BC_NODE_IS_CORE
qs = Island.add_zeros("salinity_flux", at="link")

salinity = Island.add_zeros('salinity', at='node')
Island.status_at_node[veg_array.flatten()==0] = Island.BC_NODE_IS_FIXED_VALUE
salinity[:] = 35.0
salinity[:21502]=85.000011
#salinity[:10000]=45.000011

#Island.set_status_at_node_on_edges(right=Island.BC_NODE_IS_CLOSED,
#                               top=Island.BC_NODE_IS_CLOSED,
#                               left=Island.BC_NODE_IS_CLOSED,
#                               bottom=Island.BC_NODE_IS_CLOSED)

#imshow_grid(Island, Island.status_at_node, color_for_closed='blue')


Enet = 1.15 #m/yr net evaporation rate
D = 23 #m^2/yr 
#dt = 0.005
dt = 0.2 * Island.dx * Island.dx / D
#gradients = Island.calc_grad_at_link(salinity)
#qs[Island.active_links] = -D * gradients[Island.active_links]
#dy = -Island.calc_flux_div_at_node(qs)
#salinity[Island.core_nodes] = (salinity[Island.core_nodes] +Enet)


print(dt * 100)

for i in range(10000):
    g = Island.calc_grad_at_link(salinity)
    qs[Island.active_links] = -D * g[Island.active_links]
    dqdx = Island.calc_flux_div_at_node(qs)
    dsdt = -dqdx + Enet
    salinity[Island.core_nodes] = salinity[Island.core_nodes] + (dsdt[Island.core_nodes]*dt)


test = salinity.reshape((141,305))
test1 = test[:,150]


salinityBinary = salinity
salinityBinary[salinityBinary>100]=0
salinityBinary[salinityBinary==35]=0
salinityBinary[salinityBinary==85.000011]=0
#salinityBinary[(salinityBinary > 35) | (salinityBinary < 69)] = 1
t = salinityBinary.reshape((141,305))

#Island.imshow('salinity',cmap='coolwarm')


plt.imshow(t)
plt.colorbar()
plt.contour(t,[72], colors='white')

'''
#Island.status_at_node[salinity==36]=Island.BC_NODE_IS_CORE
qs = Island.add_zeros("salinity_flux", at="link")
#Island.imshow('vegetation')
#Island.imshow('salinity')

Enet = 0.9 #m/yr net evaporation rate
D = 21 #m^2/yr 
dt = 0.2 * Island.dx * Island.dx / D

gradients = Island.calc_grad_at_link(salinity)
qs[Island.active_links] = -D * gradients[Island.active_links]
dy = -Island.calc_flux_div_at_node(qs)
salinity[Island.core_nodes] = (salinity[Island.core_nodes] + dy[Island.core_nodes]+Enet)

'''



































'''
for i in range(10):
    gradients = Island.calc_grad_at_link(salinity)
    qs[Island.active_links] = -D * g[Island.active_links]
    dzdt = -Island.calc_flux_div_at_node(qs)
    salinity[Island.core_nodes] = (salinity[Island.core_nodes] + dzdt[Island.core_nodes]+Enet) * dt

'''
#Island.imshow('salinity')

'''
for _ in range(1000):
    g = Island.calc_grad_at_link(test)
    qs[Island.active_links] = -D * g[Island.active_links]
    dzdt = -Island.calc_flux_div_at_node(qs)
    np.add(test[Island.core_nodes], ((dzdt[Island.core_nodes] - Enet) * dt),out=test[Island.core_nodes], casting="unsafe")

Island.imshow('salinity')
'''

'''
for _ in range(25):
    g = Island.calc_grad_at_link(z)
    qs[Island.active_links] = -D * g[Island.active_links]
    dzdt = -Island.calc_flux_div_at_node(qs)
    z[Island.core_nodes] += (dzdt[Island.core_nodes] - Enet) * dt
 
Island.imshow(z)



Island.imshow("salinity")

imshow_grid(Island,'salinity',at='node')

'''









'''

####turns in shapefile into a 1m raster
from osgeo import gdal

shapefile = '/Users/isamarcortes/Dropbox/Isamar/Papers_In_Prep/Paper_2/NASA_Project_Year1/G-LiHT/Regions/Polylines_All_Regions/PR1_New.shp'

# Define NoData value of new raster
NoData_value = -9999

# Filename of input OGR file
vector_fn = shapefile

# Filename of the raster Tiff that will be created
raster_fn = '/Users/isamarcortes/Dropbox/Isamar/Papers_In_Prep/Paper_4/RasterLayersForLandlab/PR1.tif'

# Open the data source and read in the extent
source_ds = gdal.OpenEx(vector_fn)
pixel_size = 0.00001  # about 25 metres(ish) use 0.001 if you want roughly 100m 

gdal.Rasterize(raster_fn, source_ds, format='GTIFF', outputType=gdal.GDT_Byte, creationOptions=["COMPRESS=DEFLATE"], noData=NoData_value, initValues=NoData_value, xRes=pixel_size, yRes=-pixel_size, allTouched=True, burnValues=1)
'''




















