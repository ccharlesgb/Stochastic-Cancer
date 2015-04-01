# -*- coding: utf-8 -*#

import math

import numpy as np
import matplotlib.pyplot as plt
import MatTools
import DeterministicTauLeap
import TauLeap
import TauLeapParam
    
def IsFixed(param):
    exists_bound = (param.r[1] > param.r[2]/(1.0-param.u[2])) and (param.r[2] < param.r[2]/(1.0-param.u[2]))
    exists_react = ((1.0-param.u[1])*param.r[0] > (1.0-param.u[2])*param.r[1]) and ((1.0-param.u[1])*param.r[0] > param.r[2])
    return exists_bound or exists_react    
    
minr2 = 0.5
maxr2 = 1.5

minr1 = 0.5
maxr1 = 1.5

mapSize = 20

cellTypes = 3
myParam = TauLeapParam.Params(cellTypes)

#Create the simulator
myGillespie = TauLeap.Sim(cellTypes)
myGillespie.params = myParam
myParam.Hook(myGillespie)
#Set History
myHist = TauLeapParam.Hist(cellTypes)
myGillespie.history = myHist

timeLimit = 1e3

myGillespie.timeLimit = timeLimit
myGillespie.RECORD_TAU_INFO = 0

myGillespie.n_c = 10
myGillespie.stopAtAppear = 0
myGillespie.epsilon = 0.05

#Param Stuff
myParam.n0[0] = 1e4

myParam.USE_D = False

myParam.r[0] = 1.0
myParam.r[1] = 1.0
myParam.r[2] = 1.0

myParam.u[0] = 0.1
myParam.u[1] = 0.1
myParam.u[2] = 0.1

#Deter Stuff
deterSim = DeterministicTauLeap.Sim(cellTypes)
deterSim.params = myParam
myHist2 = TauLeapParam.Hist(cellTypes)
deterSim.history = myHist2

deterSim.timeStep = 0.01
deterSim.timeLimit = timeLimit
deterSim.stopAtAppear = 0

SDP = 60
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
        
        myParam.r[1] = r1
        myParam.r[2] = r2
        
        if IsFixed(myParam) == 0:
            totalFixTime = 0.0
            res = myGillespie.SimulateBatch(SDP)
            avgFixTime = res.avgFixTime
            deterSim.Integrate()
            theoryTot = deterSim.curTime
        else:
            avgFixTime = float('nan')
            theoryTot = float('nan')
        
        #print("PREDICTED {0}".format(deterSim.curT))
        
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

file_name = "DetermTauLeapFracError_Size_{0}_SDP_{1}_N_{2}".format(mapSize, SDP,myParam.n0[0])

data = dict()
data["xCoords"] = x
data["yCoords"] = y
data["fracError"] = fixTimeError

MatTools.SaveDict(file_name, data)