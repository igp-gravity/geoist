#!/usr/bin/env python
# coding: utf-8

# In[239]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
from statsmodels.tsa.stattools import adfuller
from statsmodels.tsa.statespace.sarimax import SARIMAX
from statsmodels.tsa.seasonal import seasonal_decompose
from scipy.signal import detrend


# In[345]:


def despike(dt,upper_bound=None,lower_bound=None,w_size=5,*args,**kwargs):
    diff_dt = dt.diff()
    if upper_bound is None:
        upper_bound = diff_dt.std()*4
    if lower_bound is None:
        lower_bound = upper_bound / 2
    sum_diff = diff_dt.rolling(window=w_size).sum()
    is_spike = (diff_dt.abs() > upper_bound) & (sum_diff.abs() < lower_bound)
    dt[is_spike] = np.nan
    return dt.interpolate()
def dejump(dt,*args,**kwargs):
    tmp = despike(dt.diff(),*args,**kwargs)
    tmp.iloc[0] = dt.iloc[0]
    return tmp.cumsum()
def despike_v2(dt,th=None,*args,**kwargs):
    diff_dt = dt.diff()
    if th is None:
        th = diff_dt.std()*2
    is_spike = (diff_dt-diff_dt.mean()).abs() > th
    diff_dt[is_spike] = np.nan
    method = kwargs.get('method','linear')
    order = kwargs.get('order',1)
    diff_dt = diff_dt.interpolate(method=method,order=order)
    diff_dt.iloc[0] = dt.iloc[0]
    return diff_dt.cumsum()

