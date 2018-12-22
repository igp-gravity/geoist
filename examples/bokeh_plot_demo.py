# -*- coding: utf-8 -*-
"""
Created on Wed Dec 19 15:58:10 2018

@author: chens
"""

import bokeh.plotting as bp
from bokeh.palettes import Spectral4
from bokeh.models import ColumnDataSource, DatetimeTickFormatter
import numpy as np
import pandas as pd
from math import pi
#from bokeh.charts import TimeSeries
#from bokeh.io import output_file, show

d=pd.read_csv(r"D:\MyWorks\geoist\examples\data\tsdata.csv",parse_dates=True,index_col=[0])
p = bp.figure(title="lines for demo", plot_width=600, plot_height=400, x_axis_type="datetime")
#data = dict(VAL=d["value"])
#p = TimeSeries(data,y=['VAL'], title="Python Interpreters", ylabel='Stock Prices'
#               , plot_width=600, plot_height= 400, legend=True)
bp.output_file("./demo.html")
source = ColumnDataSource(d)

for data, name, color in zip(['value', 'val2', 'val3'], ["AAPL", "IBM", "MSFT"], Spectral4):
  p.line('date', data, source=source,color=color, alpha=0.8, legend=name)

p.xaxis.axis_label = 'Date'
p.yaxis.axis_label = 'Value'
p.legend.location = "top_left"
p.legend.click_policy="hide"
#p.xaxis.major_label_orientation = pi/4
bp.show(p)
#ax = np.arange(100)
#ay = np.random.rand(100)
#p = bp.figure(title="line", plot_width=600, plot_height=400)
#p.line(ax, ay)
#bp.save(p)

