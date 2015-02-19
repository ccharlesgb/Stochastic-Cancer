# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 17:39:03 2015

@author: Connor
"""

import math
import matplotlib.pyplot as plt

dataX = []
dataY = []

for i in range(0, 100):
    dataX.append(i - 50)
    dataY.append(math.exp(-float(dataX[i] * dataX[i])/100.0))
    
plt.plot(dataX, dataY)
plt.yscale("log")