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

mapSize = 30

mySim = SimTools.Gillespie(100)
mySim.r0 = 1.0
mySim.r1 = 1.0
mySim.r2 = 1.0
mySim.u1 = 0.1
mySim.u2 = 0.1
mySim.timeLimit = 100.0

simsPerDataPoint = 20

avgGoneTime = np.zeros((mapSize, mapSize))

def n0IsGone(sim):
    if sim.n0 == 0:
        sim.n0Extinct = sim.curTime
        sim.curTime = sim.timeLimit
    else:
        sim.n0Extinct = sim.timeLimit

mySim.preSim = n0IsGone

for ir1 in range(0, mapSize):
    for ir2 in range(0, mapSize):
        r1 = (float(ir1)/(mapSize-1))*(maxr1-minr1)+minr1
        r2 = (float(ir2)/(mapSize-1))*(maxr2-minr2)+minr2
        
        print("r1 = {0} r2 = {1}".format(r1,r2))     
        
        mySim.r1 = r1
        mySim.r2 = r2
        
        totalGoneTime = 0.0
        for sim in range(0,simsPerDataPoint):
            mySim.Simulate()
            if mySim.n0 == 0:
                totalGoneTime += mySim.n0Extinct
            else:
                totalGoneTime += mySim.timeLimit
        
        avgGoneTime[ir2, ir1] = totalGoneTime / float(simsPerDataPoint)

dy = (maxr2 - minr2) / mapSize
dx = (maxr1 - minr1) / mapSize

y, x = np.mgrid[slice(minr2, maxr2 + dy, dy),
slice(minr1, maxr1 + dx, dx)]

plt.pcolormesh(x,y,avgGoneTime)
plt.colorbar()

plt.xlabel("r1")
plt.ylabel("r2")

theoryX = []
theoryDiag = []
theoryXReac = []
theoryReac = []
theorySaddle = []

for ir1 in range(0,mapSize):
    r1 = float(ir1)/(mapSize - 1.0)*(maxr1-minr1)+minr1
    r2 = (1-mySim.u2)*r1
    theoryX.append(r1)
    theoryDiag.append(r2)

r1_reac = mySim.r0 * (1.0 - mySim.u1)/(1.0 - mySim.u2)
theoryXReac.append(r1_reac)
theoryXReac.append(r1_reac)
theoryReac.append(minr2)
theoryReac.append((1.0 - mySim.u2) * r1_reac)

theoryXReac.append(np.nan)
theoryReac.append(np.nan)

r2_reac = (1.0 - mySim.u1) * mySim.r0   
theoryXReac.append(minr1)
theoryXReac.append(r2_reac / (1.0 - mySim.u2))
theoryReac.append(r2_reac)
theoryReac.append(r2_reac)


plt.plot(theoryX, theoryDiag, 'k')
plt.plot(theoryXReac, theoryReac, 'k')
plt.show()

plt.xlim(minr1,maxr1)
plt.ylim(minr2,maxr2)