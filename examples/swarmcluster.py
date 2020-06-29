# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 11:10:46 2020

@author: chens
"""

from sklearn.cluster import DBSCAN
import numpy as np
X = np.array([[1, 2], [2, 2], [2, 3],
...               [8, 7], [8, 8], [25, 80]])
clustering = DBSCAN(eps=3, min_samples=2).fit(X)
clustering.labels_
clustering