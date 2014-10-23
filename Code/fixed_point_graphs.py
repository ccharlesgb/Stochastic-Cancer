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
    
minr1 = 0.5
maxr1 = 1.5

minr2 = 0.5
maxr2 = 1.5

mapSize = 50

mySim = SimTools.Gillespie(100)
mySim.r0 = 1.0
mySim.r1 = 1.0
mySim.r2 = 1.0
mySim.u1 = 0.1
mySim.u2 = 0.1
mySim.timeLimit = 100.0

simsPerDataPoint = 50

avgFixTime = np.zeros((mapSize, mapSize))

xticks = np.arange(0, mapSize, 5.0/(mapSize-1.0))
xlabels = np.arange(minr1, maxr1, 1.0/(mapSize-1.0))

yticks = np.arange(0, mapSize, 5.0/(mapSize-1.0))
ylabels = np.arange(minr2, maxr2, 1.0/(mapSize-1.0))

for ir1 in range(0, mapSize):
    for ir2 in range(0, mapSize):
        r1 = (float(ir1)/(mapSize-1))*(maxr1-minr1)+minr1
        r2 = (float(ir2)/(mapSize-1))*(maxr2-minr2)+minr2
        
        print("r1 = {0} r2 = {1}".format(r1,r2))        
        
        mySim.r1 = r1
        mySim.r2 = r2
        
        totalFixTime = 0.0
        for sim in range(0,simsPerDataPoint):
            mySim.Simulate()
            totalFixTime += mySim.curTime
        
        avgFixTime[ir2, ir1] = totalFixTime / float(simsPerDataPoint)

plt.pcolormesh(avgFixTime)
plt.colorbar()
plt.xticks(xticks, xlabels)
plt.yticks(yticks, ylabels)
plt.xlabel("r1")
plt.ylabel("r2")

theoryX = []
theoryDiag = []
theoryReac = []
theorySaddle = []

for ir1 in range(0,mapSize):
    r1 = (float(ir1))*(maxr1-minr1)+minr1
    r2 = (1-mySim.u2)*r1
    theoryX.append(r1)
    theoryDiag.append(r2)
    
plt.plot(theoryX, theoryDiag, 'k')
plt.show()