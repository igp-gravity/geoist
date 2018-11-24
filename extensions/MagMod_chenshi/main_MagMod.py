"""
    此代码分为四个部分，第一部分为参数替换区，本部分代码的目的为接入节点页面中设置的参数，
并将之替换为代码中使用的参数变量。用户在使用时，可以只修改此处的参数，不用去修改代码中的参数。
    第二部分为数据转换区，功能为接入前节点传入的数据，并转化为代码可以处理的格式。
	第三部分为代码执行区，用户可以在此编辑将要执行的代码。在实际使用时可以重点关注此部分。
	第四部分为结果导出区，此部分的功能为导出符合之后节点使用的结果文件。

"""

import pandas as pd
from pandas import DataFrame
import numpy as np                                                  
from shutil import copyfile

#接入节点界面中所设参数，并转换成固定代码。
magmod_str = '[$magmod$]'              #地磁模型IGRF11;IGRF12,分别赋值0;1,字符型
unit_str = '[$unit$]'                  #数据单位nT,赋值0，字符型
magmod = int(magmod_str)                 #转换数据类型
unit = int(unit_str)                     #转换数据类型


#获取前节点数据
basepath = GetNodeBasePath()    #获取自定义节点路径
df = CsvFile2DataFrame(basepath + 'magdata.txt')   

lon=df['Longitude']
lat=df['Latitude']
elev=df['Elevation']

#处理逻辑接口-仅供测试
val=[]
for i in range(0,len(lon)):
   if magmod == 0:
     xx=lon[i]*2+lat[i]*0.1+elev[i]*0.8
   else:
     xx=lon[i]*0.2+lat[i]*0.1+elev[i]*0.8  
   val.append(xx)

#输出接口
df1=pd.DataFrame({"A":lon,
                  "B":lat,
                  "C":elev,
                  "D":val})
OutputTable(df1, basepath + 'MagModOutPut.txt')
