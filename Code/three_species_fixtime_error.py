# -*- coding: utf-8 -*#

import TwoSpecies
import math
import anfixtime

import SimTools

import numpy as np
import matplotlib.pyplot as plt
    
minr2 = 0.5
maxr2 = 1.5

minu2 = 0.01
maxu2 = 0.2

mapSize = 10

sim3 = SimTools.Gillespie(10)
sim3.timeLimit = 100000
sim3.u1 = 0.1
sim3.u2 = 0.01
sim3.in0 = 10

sim3.r0 = 1.5
sim3.r1 = 1.0
sim3.r2 = 1.0

#This simulator considers the transition from type0 to type1
sim01 = TwoSpecies.Gillespie(10)
sim01.u1 = sim3.u1

sim01.r0 = sim3.r0
sim01.r1 = sim3.r1

sim01.ij = 0

#Now consider transition from type1 to type2
sim12 = TwoSpecies.Gillespie(10)
sim12.u1 = sim3.u2

sim12.r0 = sim3.r1
sim12.r1 = sim3.r2

sim01.ij = 0

simsPerDataPoint = 200

fixTimeError = np.zeros((mapSize, mapSize))

xticks = np.arange(0, mapSize, 5.0/(mapSize-1.0))
xlabels = np.arange(minu2, maxu2, 1.0/(mapSize-1.0))

yticks = np.arange(0, mapSize, 5.0/(mapSize-1.0))
ylabels = np.arange(minr2, maxr2, 1.0/(mapSize-1.0))

for ir2 in range(0, mapSize):
    for iu2 in range(0, mapSize):
        r2 = (float(ir2)/(mapSize-1))*(maxr2-minr2) + minr2
        u2 = (float(iu2)/(mapSize-1))*(maxu2-minu2) + minu2
        
        print("r2 = {0} u2 = {1}".format(r2,u2))        
        
        sim3.u2 = u2
        sim3.r2 = r2
        
        totalFixTime = 0.0
        for sim in range(0,simsPerDataPoint):
            sim3.Simulate()
            totalFixTime += sim3.curTime
            if sim3.Fixated() == 0:
                print("WARNING u2 = {0} r2 = {1} NO FIXATE".format(u2,r2))
        
        sim01.u1 = sim3.u1
        sim12.u1 = sim3.u2
        
        sim12.r1 = sim3.r2
    
        theory01 = anfixtime.GetFixTimeJ(sim01, sim01.ij)
        theory12 = anfixtime.GetFixTimeJ(sim12, sim12.ij)
    
        theoryTot = theory01 + theory12   

        avgFixTime = (totalFixTime / float(simsPerDataPoint))        
        
        fixTimeError[iu2, ir2] = (theoryTot - avgFixTime) / avgFixTime
        print("FixTimeErr {0}".format(fixTimeError[iu2, ir2]))

        
dy = (maxr2 - minr2) / mapSize
dx = (maxu2 - minu2) / mapSize

y, x = np.mgrid[slice(minr2, maxr2 + dy, dy),
slice(minu2, maxu2 + dx, dx)]

plt.pcolormesh(x,y,fixTimeError)
plt.colorbar()

#plt.xticks(xticks, xlabels)
#plt.yticks(yticks, ylabels)

plt.xlabel("u2")
plt.ylabel("r2")


plt.xlim(minu2,maxu2)
plt.ylim(minr2,maxr2)