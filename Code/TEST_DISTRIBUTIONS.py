# -*- coding: utf-8 -*-
"""
Created on Sat May 02 16:02:10 2015

@author: Connor
"""

import wright_fisher
import matplotlib.pyplot as plt
import math
import SimUtil
import TauSolver
import TauLeapParam
import MatTools
import random
import numpy

def LogNormalU2(mean,variance, count):
    log_mean = math.log10(mean)
    u = []
    for i in range(0,count):
        log_u = random.normalvariate(log_mean, variance)
        u.append(math.pow(10.0, log_u))
    return u
    
data = []
COUNT = int(1e3)

meanU = 1e-6
variance = 0.2
cellTypes = 21

for i in range(0, COUNT):
    u_arr = LogNormalU2(meanU, variance, cellTypes-1)
    
    for i in range(0, cellTypes-1):
        data.append(math.log10(u_arr[i]))
        
print("MYPARAM U ", min(data), max(data))

plt.hist(data,bins=100)