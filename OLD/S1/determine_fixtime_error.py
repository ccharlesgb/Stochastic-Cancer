# -*- coding: utf-8 -*#

import TwoSpecies
import math
import anfixtime

import SimTools

import numpy as np
import matplotlib.pyplot as plt

import MatTools

import DetermenisticModel
    
minr2 = 0.5
maxr2 = 1.5

minr1 = 0.5
maxr1 = 1.5

mapSize = 8

sim3 = SimTools.Gillespie(10)
sim3.timeLimit = 100000
sim3.u1 = 0.1
sim3.u2 = 0.1
sim3.in0 = 10

sim3.r0 = 1.0
sim3.r1 = 1.0
sim3.r2 = 1.0

deterSim = DetermenisticModel.Determenistic()

simsPerDataPoint = 1
fixTimeError = np.zeros((mapSize, mapSize))

xticks = np.arange(0, mapSize, 5.0/(mapSize-1.0))
xlabels = np.arange(minr1, maxr1, 1.0/(mapSize-1.0))

yticks = np.arange(0, mapSize, 5.0/(mapSize-1.0))
ylabels = np.arange(minr2, maxr2, 1.0/(mapSize-1.0))

for ir1 in range(0, mapSize):
    for ir2 in range(0, mapSize):
        r1 = (float(ir1)/(mapSize-1))*(maxr1-minr1)+minr1
        r2 = (float(ir2)/(mapSize-1))*(maxr2-minr2)+minr2
        #u2 = (float(iu2)/(mapSize-1))*(maxu2-minu2) + minu2
        
        print("r1 = {0} r2 = {1}".format(r1,r2))        
        
        sim3.r1 = r1
        sim3.r2 = r2
        
        if deterSim.IsFixed(sim3) == 0:       
            totalFixTime = 0.0
            for sim in range(0,simsPerDataPoint):
                sim3.Simulate()
                totalFixTime += sim3.curTime
                if sim3.Fixated() == 0:
                    print("WARNING r1 = {0} r2 = {1} NO FIXATE".format(r1,r2))
        else:
            totalFixTime = float('nan')
        
        deterSim.Reset()
        deterSim.Integrate(sim3)
        #print("PREDICTED {0}".format(deterSim.curT))
        theoryTot = deterSim.curT

        avgFixTime = (totalFixTime / float(simsPerDataPoint))        
        
        fixTimeError[ir2, ir1] = (theoryTot - avgFixTime) / avgFixTime
        print("FixTimeErr {0}".format(fixTimeError[ir1, ir2]))

        
dy = (maxr2 - minr2) / mapSize
dx = (maxr1 - minr1) / mapSize

y, x = np.mgrid[slice(minr2, maxr2 + dy, dy),
slice(minr1, maxr1 + dx, dx)]

plt.pcolormesh(x,y,fixTimeError)
plt.colorbar()

#plt.xticks(xticks, xlabels)
#plt.yticks(yticks, ylabels)

plt.xlabel("r1")
plt.ylabel("r2")

plt.clim([-1.0,1.0])

plt.xlim(minr1,maxr1)
plt.ylim(minr2,maxr2)

file_name = "three_species_errordeterm_Size_{0}_SDP_{1}_N_{2}".format(mapSize, simsPerDataPoint,sim3.in0)
MatTools.ColourMap(file_name, x, y, fixTimeError, xLabel = "r1", yLabel = "r2", zLabel = "frac_error")