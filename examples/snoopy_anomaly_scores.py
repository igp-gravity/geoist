# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 19:09:03 2018
异常检测模块
Calculate anomaly scores
@author: chens
"""

from geoist.snoopy.anomaly_detector import AnomalyDetector

ts = {0: 0, 1: 0.5, 2: 1, 3: 1, 4: 1, 5: 0, 6: 0, 7: 0, 8: 0}

my_detector = AnomalyDetector(ts)
score = my_detector.get_all_scores()
for timestamp, value in score.iteritems():
    print(timestamp, value)