# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 22:42:29 2019

@author: chens
"""
import pathlib
import os
import requests
from geoist.others.utils import file_hash,make_registry

tmppath = pathlib.Path(__file__).parent
path = pathlib.Path(tmppath,'data')
file1 = pathlib.Path(path,'cncat_mw.csv')
file2 = pathlib.Path(path,'isc-gem-cat.csv')
print(file_hash(file1))
print(file_hash(file2))
regfile = pathlib.Path(path,'regHash.txt')

make_registry(path, regfile, recursive=True)

import pandas as pd

dd=pd.read_csv(r'D:\MyWorks\geoist\examples\data\ttt\cncat_mw.txt', sep='\t', encoding = 'utf-8')
##先获得HASH信息从Github到本地
#url = "https://raw.githubusercontent.com/gravity-igpcea/dataset/master/"
#filename = url + 'regHash.txt'
#downfile = pathlib.Path(path,'regHashd.txt')
#f = open(downfile,'wb+')
#with f:
#    r = requests.get(filename, stream=True)
#    r.raise_for_status()
#    for chunk in r.iter_content(chunk_size=1024):
#        if chunk:
#            f.write(chunk)
##print(r.content)
##f.write(str(r.content)
#from geoist.others import fetch_data
#fetch_data.Drepo1.load_registry(downfile)
##下载数据
#print(os.path.join(os.path.dirname(__file__)))
#cncat = fetch_data.fetch_catalog()  #dataframe