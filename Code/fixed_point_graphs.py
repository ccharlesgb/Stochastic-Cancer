# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np

import SimTools

def FixedPointExists(r1, r2, u2):
    condition = r2 / (1.0 - u2)
    return r1 > condition and r2 < condition
    
def FixedPointStable(r0, r1, u2):
    return (1.0 - u2)*r1 > (1.0 - u2) * r0
    
def FixedPointSaddle(r0,r1, u1, u2):
    return (1.0 - u2)*r1 < (1.0 - u1)*r0
    
minr1 = 0.0
maxr1 = 2.0

minr2 = 0.0
maxr2 = 2.0

mapSize = 10

mySim = SimTools.Gillespie(10)
mySim.r0 = 1.0
mySim.r1 = 1.0
mySim.r2 = 1.0
mySim.u1 = 0.1
mySim.u2 = 0.1
mySim.timeLimit = 100000.0

simsPerDataPoint = 10

avgFixTime = np.zeros((mapSize, mapSize))

for ir1 in range(0, mapSize):
    for ir2 in range(0, mapSize):
        r1 = (float(ir1)/(mapSize-1))*(maxr1-minr1)+minr1
        r2 = (float(ir2)/(mapSize-1))*(maxr2-minr2)+minr2
        
        mySim.r1 = r1
        mySim.r2 = r2
        
        totalFixTime = 0.0
        for sim in range(0,simsPerDataPoint):
            mySim.Simulate()
            totalFixTime += mySim.curTime
        
        avgFixTime[ir1, ir2] = totalFixTime / float(simsPerDataPoint)

plt.pcolor(avgFixTime)