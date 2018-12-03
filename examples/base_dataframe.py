# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 09:24:41 2018
测试datafram的操作
@author: chens
"""

import pandas as pd
import numpy as np
#读入数据接口
df1=pd.read_csv("D:\MyWorks\geoist-master\pfm\gradata.csv")
lon=df1['Longitude']
lat=df1['Latitude']
elev=df1['Elevation']
tidex=0 #or 1 2 
unit=0 #or 1 单位类型
newindex=0 #or 1 是否新建列
#处理逻辑接口-仅供测试
val=[];
for i in range(0,len(lon)):
   if unit == 0:
     xx=lon[i]*2+lat[i]*0.1+elev[i]*0.8
   else:
     xx=lon[i]*0.2+lat[i]*0.1+elev[i]*0.8  
   val.append(xx)
   
#输出接口
export={}
export['lon']=lon
export['lat']=lat
export['elev']=elev
export['val']=val
df2=pd.DataFrame(export)
df2.to_csv("D:\MyWorks\geoist-master\pfm\df2.csv")