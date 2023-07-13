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

'''
dx = 1 ###space per cell
rows = 40###rows of grid
columns = 40###columns of grid

grid_veg = RasterModelGrid((rows,columns),dx)###creating grid
veg = grid_veg.add_ones('vegetation',at='cell')###initializing vegetation grid




####salinity grid

grid_salinity = RasterModelGrid((rows,columns),dx)###creating same sized grid for salinity
salinity = grid_salinity.add_zeros('salinity',at='node',units = 'ppt',clobber = True)###initializing grid with zeros

salinityCenterX,salinityCenterY = (20,20)###position of highest outer edge salinity
salinity_node = int(salinityCenterY / dx * columns + salinityCenterX / dx)###this part came from the tutorial

grid_salinity.at_node['salinity'][salinity_node]=130###high in center of island
Enet = 2 ###source term for model
D=100
dx, dy = grid_salinity.dx, grid_salinity.dy ###getting dx and dy from rasterModelGrid
dt = 1000
n_steps = 100 ###total number of iterations
salinityAtCenter = np.zeros(n_steps)
'''
#imshow_grid(grid_salinity,'salinity')

#ld = LinearDiffuser(grid_salinity,linear_diffusivity=10)


#Attempting to work with NetworkModelGrid 
#Past this point

import rasterio as rio





File = rio.open('/Users/isamarcortes/Dropbox/Isamar/Papers_In_Prep/Paper_4/RasterLayersForLandlab/PR1.tif')
veg_array = File.read(1)
sal_array = File.read(1)
sal_array[sal_array==0]=35
sal_array[sal_array==1]=36



#plt.imshow(array)
#plt.colorbar()

Island = RasterModelGrid((veg_array.shape))
Island.add_field('vegetation',veg_array,at='node')



test = Island.add_field('salinity',sal_array,at='node')
Island.status_at_node['salinity'==35]=Island.BC_NODE_IS_CLOSED
Island.status_at_node['salinity'==36]=Island.BC_NODE_IS_CORE
qs = Island.add_zeros("salinity_flux", at="link")
#Island.imshow('vegetation')

Enet = 0.9 #m/yr net evaporation rate
D = 21 #m^2/yr 
dt = 0.2 * Island.dx * Island.dx / D


#Island.set_closed_boundaries_at_grid_edges(True, False, True, False)

#g = Island.calc_grad_at_link(test)
#q = -D * g
#dqda = Island.calc_flux_div_at_node(q) 

for i in range(10):
    g = Island.calc_grad_at_link(test)
    qs[Island.active_links] = -D * g[Island.active_links]
    dzdt = -Island.calc_flux_div_at_node(qs)
    test[Island.core_nodes] = (test[Island.core_nodes]+ dzdt[Island.core_nodes]+Enet) * dt


Island.imshow('salinity')

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




















