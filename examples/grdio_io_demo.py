
# coding: utf-8

# In[1]:

import numpy as np
from geoist import pfm


# In[2]:

# load original dataset
grd1 = pfm.grdio.grddata()
grd1.load_surfer('.\\data\\Demogrid.grd')


# In[3]:

# test export binary/ascii
grd1.export_surfer('.\\data\\Demogrid_b.grd',file_format='b')
grd1.export_surfer('.\\data\\Demogrid_a.grd',file_format='a')


# In[4]:

# test load binary surfer
grd2 = pfm.grdio.grddata()
grd2.load_grd('.\\data\\Demogrid_b.grd')


# In[5]:

# test load ascii surfer
grd3 = pfm.grdio.grddata()
grd3.load_grd('.\\data\\Demogrid_a.grd')


# In[6]:

# test load ascii grd
grd4 = pfm.grdio.grddata()
grd4.load_grd('.\\data\\GridFile1.grd')

