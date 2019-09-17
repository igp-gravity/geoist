
# coding: utf-8

# In[1]:

import pandas as pd
import matplotlib.pyplot as plt
get_ipython().magic('matplotlib inline')


# In[2]:

# load data
line_data = pd.read_csv('line_data.csv',index_col=0,parse_dates=True)
annotation = pd.read_csv('annotation.csv',index_col=0,parse_dates=True)


# In[3]:

# get annotation's coordinate
ind = line_data.index.union(annotation.index)
d2 = line_data.reindex(ind).interpolate(method='linear').fillna(0)


# In[4]:

# set upper and lower threshhold
upper = 0.75
lower =-0.75


# In[5]:

# draw the line_data
ax = line_data.plot()

# fill between line_data and threshhold
ax.fill_between(x=line_data.index,y1=line_data.iloc[:,0],y2=upper,where=line_data.iloc[:,0]>upper,color='salmon')
ax.fill_between(x=line_data.index,y1=line_data.iloc[:,0],y2=lower,where=line_data.iloc[:,0]<lower,color='skyblue')

# set annotation properties
sep = line_data.max()-line_data.min()
arrowprops = {'arrowstyle':'->'}
# draw annotation
for row in annotation.itertuples():
    ax.annotate(row[1],xy=(row[0],d2.loc[row[0]].iloc[0]),xytext=(row[0],d2.loc[row[0]].iloc[0]+0.1*sep),annotation_clip=True,ha='center',arrowprops=arrowprops)

