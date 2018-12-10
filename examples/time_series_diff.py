# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 19:20:11 2018

@author: chens
"""

import pandas as pd
time_series = pd.Series([2,4,3,5,6,7,4,5,6,3,2,4], 
                        index=pd.date_range(start='2000', periods=12, freq='a'))
time_series_diff = time_series.diff(1).dropna() #差分结果

time_series_restored = pd.Series([time_series[0]], 
                                 index=[time_series.index[0]]).append(time_series_diff).cumsum()

print(time_series_diff, time_series_restored)

