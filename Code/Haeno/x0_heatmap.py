# -*- coding: utf-8 -*-
"""
Created on Sun Nov 09 15:52:45 2014

@author: Connor
"""

#Haeno Analysis Testing
#Sweeps r1

import math
import matplotlib.pyplot as plt
import HaenoModel
import numpy as np
#import ./SimTools

myHaeno = HaenoModel.HaenoModel()

myHaeno.N = 100
myHaeno.t = 5.0
myHaeno.r0 = 1.0
myHaeno.r1 = 1.0
myHaeno.r2 = 1.0

myHaeno.u1 = 0.1
myHaeno.u2 = 0.01

myHaeno.clampB = 1

myHaeno.V_accuracy = 1e-4
myHaeno.V_maxIter = 500

mapSize = 20

minr1 = 0.5
maxr1 = 1.5

minr2 = 0.5
maxr2 = 1.5

x0Prob = np.zeros((mapSize, mapSize))

xticks = np.arange(0, mapSize, 5.0/(mapSize-1.0))
xlabels = np.arange(minr1, maxr1, 1.0/(mapSize-1.0))

yticks = np.arange(0, mapSize, 5.0/(mapSize-1.0))
ylabels = np.arange(minr2, maxr2, 1.0/(mapSize-1.0))

for ir1 in range(0, mapSize):
    for ir2 in range(0, mapSize):
        r1 = (float(ir1)/(mapSize-1))*(maxr1-minr1)+minr1
        r2 = (float(ir2)/(mapSize-1))*(maxr2-minr2)+minr2
        
        print("r1 = {0} r2 = {1}".format(r1,r2)) 
        
        myHaeno.r1 = r1
        myHaeno.r2 = r2
        
        x0Prob[ir2, ir1] = myHaeno.GetX0()

dy = (maxr2 - minr2) / mapSize
dx = (maxr1 - minr1) / mapSize

y, x = np.mgrid[slice(minr2, maxr2 + dy, dy),
slice(minr1, maxr1 + dx, dx)]

plt.pcolormesh(x,y,x0Prob)
plt.colorbar()

plt.xlabel("r1")
plt.ylabel("r2")