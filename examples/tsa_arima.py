# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 12:29:08 2020

@author: chens
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import geoist.snoopy.tsa as tsa

# parameters for loading data
data_path = Path(tsa.__file__).parent
orig_file = Path(data_path,"data",'50002_1_2312.txt')
water_file = Path(data_path,"data",'water_level_res.txt')


# parameters for processing data
na_values = None


# load data
data=pd.read_csv(Path(orig_file),parse_dates=[[0,1]],header=None,delim_whitespace=True,index_col=[0],na_values=na_values)
data.index.name = 'time'
data.columns = ['origin_data']
ax=data.plot(figsize=(15,12),y=data.columns[0])
ax.set_xlabel('Date')
ax.set_ylabel('Value')

# despike
thresh_hold = 200.0
data['despiked'],data['flag'] = tsa.despike_v2(data['origin_data'].interpolate(),th=thresh_hold)
ax=data.plot(figsize=(15,12),y=data.columns[:-1])
ax.set_xlabel('Date')
ax.set_ylabel('Value')
plt.grid()
plt.legend()
plt.title("The preliminary result by threshold={}".format(thresh_hold),loc='left')

# detrend
data['detrend'] = tsa.detrend(data['despiked'])
data.plot(figsize=(15,12),y=['detrend','despiked'])

# ARIMA
p = 5
d = 0
q = 1
P = 0
D = 0
Q = 0
s = 0
model = tsa.SARIMAX(data['detrend'].dropna(),
                    order=(p,d,q),
                    seasonal_order=(P,D,Q,s),
                    enforce_stationarity=False)
rests = model.fit()

pred = rests.get_forecast(180)
pci = pred.conf_int()
fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, ncols=1,figsize=(12,12))
ax0.plot(data['detrend'].dropna())
ax1.plot(data['detrend'].dropna().values)
pred.predicted_mean.plot(ax=ax1,label='forecast')
ax1.fill_between(pci.index,pci.iloc[:,0],pci.iloc[:,1],color='k',alpha=0.2,label='0.05 confidence interval')
ax1.legend()

start_day=365*10
pred = rests.get_prediction(start_day)
pci = pred.conf_int()
pred.predicted_mean.plot(ax=ax2,label='prediction')
ax2.fill_between(pci.index,pci.iloc[:,0],pci.iloc[:,1],color='k',alpha=0.2,label='0.05 confidence interval')
data['detrend'].dropna().iloc[start_day:].asfreq('5d').plot(style='o',ax=ax2)
ax2.legend()


# date_format = None

# seasonal decomposition
window_size = 50
period = 365
na_values_output = np.nan

decomposition = tsa.seasonal_decompose(data['despiked'],freq=period,extrapolate_trend='freq')

fig=decomposition.plot()
fig.set_size_inches(15,8)

data['trend'] = decomposition.trend.fillna(na_values_output)
data['seasonal'] = decomposition.seasonal.fillna(na_values_output)
data['residual'] = decomposition.resid.fillna(na_values_output)



from geoist.snoopy.anomaly_detector import AnomalyDetector
ts = dict(zip(range(len(data)),data['despiked'].values))

window_size = 10
center = False
data['ma'] = data['detrend'].rolling(window=window_size,center=center,min_periods=1).mean()

tsb = dict(zip(range(len(data)),data['ma'].values))
#my_detector = AnomalyDetector(ts, baseline_time_series = tsb, algorithm_name = 'diff_percent_threshold',
#                              algorithm_params = {'percent_threshold_upper': 20, 'percent_threshold_lower': -20}) #'exp_avg_detector'

# anomaly below baseline
algorithm_params = {'percent_threshold_upper': 20,
                    'offset': 20000,
                    'scan_window': 24,
                    'confidence': 0.01}

my_detector = AnomalyDetector(ts, baseline_time_series = tsb, algorithm_name = 'sign_test',
                               algorithm_params = algorithm_params) #'exp_avg_detector'

score = my_detector.get_all_scores()

fig = plt.figure(figsize=(15,12))
ax = fig.add_subplot(211)
ax.plot(data['despiked'].dropna().values)
ax.plot(data['ma'].dropna().values)
ax = fig.add_subplot(212)
ax.plot(score.timestamps,score.values)


