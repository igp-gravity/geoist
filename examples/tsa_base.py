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


# test stationary
window_size = 50

res = tsa.adfuller(data['despiked'].values)
tsa.print_adf(res,'despiked data')
data['mean'] = data['despiked'].rolling(window=window_size).mean()
data['std'] = data['despiked'].rolling(window=window_size).std()
data.plot(figsize=(15,12),y=['despiked','mean','std'])

# detrend
data['detrend'] = tsa.detrend(data['despiked'])
data.plot(figsize=(15,12),y=['detrend','despiked'])


# test stationary again
res = tsa.adfuller(data['detrend'].values)
tsa.print_adf(res,'detrended data')
data['mean'] = data['detrend'].rolling(window=window_size).mean()
data['std'] = data['detrend'].rolling(window=window_size).std()
data.plot(figsize=(15,12),y=['detrend','mean','std'])

# seasonal decomposition
period = 365
na_values_output = np.nan

decomposition = tsa.seasonal_decompose(data['despiked'],freq=period,extrapolate_trend='freq')

fig=decomposition.plot()
fig.set_size_inches(15,8)

data['trend'] = decomposition.trend.fillna(na_values_output)
data['seasonal'] = decomposition.seasonal.fillna(na_values_output)
data['residual'] = decomposition.resid.fillna(na_values_output)

# test stationary on residual
res = tsa.adfuller(data['residual'].dropna().values)
tsa.print_adf(res,'residual data')
data['mean'] = data['residual'].rolling(window=window_size).mean()
data['std'] = data['residual'].rolling(window=window_size).std()
data.plot(figsize=(15,12),y=['residual','mean','std'])


# first order difference
data['diff'] = data['origin_data'].diff()
ax = data.plot(figsize=(15,12),y=['origin_data','diff'])

# second order difference
data['diff2'] = data['diff'].diff()
ax = data.plot(figsize=(15,12),y=['origin_data','diff2'])

# moving average
window_size = 50
center = False
data['ma'] = data['residual'].rolling(window=window_size,center=center,min_periods=1).mean()
data.plot(y=['residual','ma'])

# exponential moving average
factor = 0.3
data['ewm'] = data['residual'].ewm(alpha=factor).mean()
data.plot(y=['residual','ewm'])

# 距平分析
# load data
dateparser = lambda x: pd.to_datetime(x,format='%Y%m')
water = pd.read_csv(water_file,header=None,parse_dates=True,index_col=0,delim_whitespace=True,date_parser=dateparser)
water[water == 99999] = np.nan
water = water.interpolate()
water.columns = ['origin','mean','departure']
water_origin = pd.DataFrame(water[water.columns[0]]).copy()

# call departure, and plot.
water_origin,_ = tsa.despike_v2(water_origin,th=200)
wate_departure = tsa.departure(water_origin)
ax = wate_departure.plot(figsize=(16,9))
ax.invert_yaxis()

# old result
# ax = water.plot(figsize=(16,9))
# ax.invert_yaxis()


# upsample
water_daily = water_origin.resample('D').asfreq().interpolate()
water_daily.head(10)

# downsample
water_monthly = water_daily.resample('MS').asfreq().interpolate()
water_monthly.head(10)

# monthly mean
water_monthly1 = water_daily.resample('MS').mean().interpolate()
water_monthly1.head(10)


# filter
# generate dataset
sample_rate = 30.0
n = np.arange(300)
orig_data = np.sin(0.1*np.pi*n)+2.0*np.cos(0.5*np.pi*n)+1.5*np.sin(0.8*np.pi*n)

# generate filter
order = 10
nyq = 0.5*sample_rate
lower_cut_rate = 7.0 / nyq
upper_cut_rate = 10.0 / nyq
sos = tsa.butter(10,lower_cut_rate,btype='low',output='sos')

# apply filter to data
filtered_data = tsa.sosfiltfilt(sos,orig_data)

# plot data
fig = plt.figure(figsize=(16,16))
ax = plt.subplot(211)
ax.plot(n/sample_rate,orig_data,label='orig_signal')
ax.plot(n/sample_rate,filtered_data,label='filtered_signal')
ax.set_xlabel('time(s)')
ax.set_ylabel('magnitude')
ax.legend()
plt.title('Effect of low pass filter (critical frequency:{}Hz)'.format(lower_cut_rate*nyq),loc='left')
plt.grid()

ax = plt.subplot(212)
w,h = tsa.sosfreqz(sos,worN=1000)
ax.plot(0.5*sample_rate*w/np.pi,np.abs(h))
ax.set_xlabel('frequency(Hz)')
ax.set_ylabel('response')
plt.title('Frequency response',loc='left')
plt.grid()

# psd
f_w,pxx_w = tsa.welch(orig_data,sample_rate,nperseg=256,scaling='spectrum')
f_p,pxx_p = tsa.periodogram(orig_data,sample_rate,scaling='spectrum')
f_l = np.linspace(0.1,14,3000)*np.pi*2.0
pxx_l = tsa.lombscargle(n/sample_rate,orig_data,f_l)

# plot result
fig = plt.figure(figsize=(15,9))
ax = fig.add_subplot(111)
ax.plot(f_w,pxx_w,label='welch')
ax.scatter(f_p,pxx_p,label='peridogram',c='g')
ax.plot(0.5*f_l/np.pi,np.sqrt(pxx_l*4.0/len(orig_data)),alpha=0.7,label='lombscargle')
ax.legend()

# ARIMA
p = 5
d = 0
q = 1
P = 0
D = 0
Q = 0
s = 0
model = tsa.SARIMAX(data['residual'].dropna(),
                    order=(p,d,q),
                    seasonal_order=(P,D,Q,s),
                    enforce_stationarity=False)
res = model.fit()


pred = res.get_forecast(180)
pci = pred.conf_int()
res111 = data['residual'].dropna().values
fig = plt.figure(figsize=(15,12))
ax = fig.add_subplot(111)
pred.predicted_mean.plot(ax=ax,label='forecast')
ax.plot(data['residual'].dropna().values)
ax.fill_between(pci.index,pci.iloc[:,0],pci.iloc[:,1],color='k',alpha=0.2,label='0.05 confidence interval')
ax.plot(res111)
ax.legend()


start_day=365*7
pred = res.get_prediction(start_day)
pci = pred.conf_int()
fig = plt.figure(figsize=(15,12))
ax = fig.add_subplot(111)
pred.predicted_mean.plot(ax=ax,label='forecast')
ax.fill_between(pci.index,pci.iloc[:,0],pci.iloc[:,1],color='k',alpha=0.2,label='0.05 confidence interval')
data['residual'].dropna().iloc[start_day:].asfreq('5d').plot(style='o',ax=ax)
ax.legend()



# date_format = None

from geoist.snoopy.anomaly_detector import AnomalyDetector
ts = dict(zip(range(len(data)),data['residual'].values))
my_detector = AnomalyDetector(ts)



score = my_detector.get_all_scores()


fig = plt.figure(figsize=(15,9))
ax = fig.add_subplot(211)
data.plot(ax=ax,y=['residual'])
ax = fig.add_subplot(212)
ax.plot(score.timestamps,score.values)


