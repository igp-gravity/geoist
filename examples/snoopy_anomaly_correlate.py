# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 19:12:25 2018
计算异常相关
Correlate ts1 with ts2 on every anomaly.
@author: chens
"""

from geoist.snoopy.anomaly_detector import AnomalyDetector
from geoist.snoopy.correlator import Correlator
import matplotlib.pyplot as plt

ts1 = {0: 0, 1: 0.5, 2: 1, 3: 1, 4: 1, 5: 0, 6: 0, 7: 0, 8: 0}
ts2 = {0: 0, 1: 0.5, 2: 1, 3: 0.5, 4: 1, 5: 0, 6: 1, 7: 1, 8: 1}

my_detector = AnomalyDetector(ts1, score_threshold=1.5)
score = my_detector.get_all_scores()
anomalies = my_detector.get_anomalies()
for a in anomalies:
    time_period = a.get_time_window()
    my_correlator = Correlator(ts1, ts2, time_period)
    if my_correlator.is_correlated(threshold=0.8):
        print("ts2 correlate with ts1 at time period (%d, %d)" % time_period)
        
for timestamp, value in score.iteritems():
    print(timestamp, value) 
       
fig = plt.figure(figsize=(15,9))
ax = fig.add_subplot(111)
ax.plot(score.timestamps,score.values)