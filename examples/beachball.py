# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 17:11:22 2020

@author: chens
"""

# from obspy.imaging.beachball import beachball
 
# mt = [0.91, -0.89, -0.02, 1.78, -1.55, 0.47] 
# beachball(mt, size=200, linewidth=2, facecolor='b')
# mt2 = [150, 87, 1]
 
# beachball(mt2, size=200, linewidth=2, facecolor='r')
 
# mt3 = [-2.39, 1.04, 1.35, 0.57, -2.94, -0.94]
 
# beachball(mt3, size=200, linewidth=2, facecolor='g')

# import geoist.catalog.QCreport as qc
# pathname,prefix = qc.create_figures_new(qc.qcinit())
#qc.generate_html(pathname,prefix,qc.to_show)

from os.path import dirname
from geoist.catalog import Catalogue as cat
from geoist.catalog import Exploration as exp
import geoist.catalog.QCreport as qc
import geoist.others.fetch_data as data
import geoist as gi
from geoist.others.fetch_data import _retrieve_file as downloadurl


## 下载地震目录
url = data.ispec_catalog_url
filename = '2020-03-25CENC-M4.7.dat'
localpath = downloadurl(url+filename, filename)
print(localpath) #文件路径

db2 = cat.Database('2020-03-25CENC-M4.7')
# header = ['Year', 'Month','Day','Hour','Minute','Second','Latitude','Longitude', 'Depth','MagType','MagSize','Log']
# db2.Import0(localpath, Header = header, Delimiter= ' ', flag = False)

db2.Import0('D:/MyWorks/geoist/geoist/catalog/data/us2010-2012.csv')
#header = ['time', 'latitude','longitude', 'depth','mag','magType','type']
#db2.Import0('D:/MyWorks/geoist/geoist/catalog/data/us2010-2012-1.csv') #, Header=header)
#db2.Size()
#db2.Events[2]
db2.Info()
exp.MagTimePlot(db2)

#print(gi.EXAMPLES_PATH, gi.DATA_PATH, gi.TEMP_PATH)
##gi.__verison__
##gi.log.info('test')
#pathname,prefix = qc.create_figures_new(qc.qcinit(), catalog_file = 'us2010-2012.csv')
pathname,prefix = qc.create_figures_new(qc.qcinit(), db2)


# import pandas as pd
# df = pd.DataFrame(columns=['newind','time','latitude','longitude', 'depth' ,'mag', 'magType'])
# dt = pd.DataFrame({'year': db2.Extract('Year'), 'month': db2.Extract('Month'),
#               'day': db2.Extract('Day'),'hour': db2.Extract('Hour'),'minute': db2.Extract('Minute'),
#               'second': db2.Extract('Second')} )
# df['newind'] = pd.to_datetime(dt)
# df['id'] = db2.Extract('Id') 
# df['time'] = [x.strftime("%Y-%m-%dT%H:%M:%S.%fZ") for x in df['newind'].tolist()]
# df['latitude'] = db2.Extract('Latitude') 
# df['longitude'] = db2.Extract('Longitude')
# df['depth'] = db2.Extract('Depth')
# df['mag'] = db2.Extract('MagSize')
# df['magType'] = db2.Extract('MagType')


#pathname,prefix = qc.create_figures_new(qc.qcinit())
#qc.generate_html(pathname,prefix,qc.to_show)


