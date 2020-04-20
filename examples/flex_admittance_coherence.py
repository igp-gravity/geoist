# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 16:55:28 2020

@author: chens
"""

import numpy as np
import pandas as pd
from pathlib import Path
from geoist.flex import TopoGrid, BougGrid, Project
from geoist import DATA_PATH

topofile = Path(DATA_PATH, 'Topo_NA.xyz')
bugfile = Path(DATA_PATH, 'Bouguer_NA.xyz')
# Read header (first line) of data set using pandas to get grid parameters
xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, nx, ny = \
pd.read_csv(topofile, sep='\t', nrows=0).columns[1:].values.astype(float)

# Change type of nx and ny from float to integers
nx = int(nx)
ny = int(ny)

# Read topography and bouguer anomaly data 
topodata = pd.read_csv(topofile, sep='\t', \
    skiprows=1, names=['x', 'y', 'z'])['z'].values.reshape(ny,nx)[::-1]
bougdata = pd.read_csv(bugfile, sep='\t', \
    skiprows=1, names=['x', 'y', 'z'])['z'].values.reshape(ny,nx)[::-1]

# Load the data as `TopoGrid` and `BougGrid` objects
topo = TopoGrid(topodata, dx, dy)
boug = BougGrid(bougdata, dx, dy)

# Create contours of coastlines
contours = topo.make_contours(0.)

# Make mask over deep basins
mask = (topo.data < -500.)

# Plot topo and boug objects with mask and contours - use alternative colormaps for topo
topo.plot(mask=mask, contours=contours, cmap='gist_earth', vmin=-1000, vmax=2250)
boug.plot(mask=mask, contours=contours, cmap='seismic', vmin=-400, vmax=400)

project = Project(grids=[topo, boug])

# Define empty Project
# project = Project()

# # Add topo
# project += topo

# # Add boug
# project.append(boug)
#project.init()
# project.wlet_admit_coh()
# project.__dict__.keys()
# project.plot_admit_coh(kindex=7, contours=contours, mask=mask)
# # Take random cell value within grid and set as attribute
# project.cell = (250, 100)
# # Plot admittance and coherence functions
# project.plot_functions()