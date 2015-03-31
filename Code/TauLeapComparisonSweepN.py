# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 16:23:16 2015

@author: Connor
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Mar 09 15:09:08 2015

@author: Connor
"""

import matplotlib.pyplot as plt
import MatTools
import numpy as np
import math
import time
import DeterministicTauLeap

import TauLeap
import TauLeapParam
import SimUtil

cellTypes = 3
myParam = TauLeapParam.Params(cellTypes)

#Create the simulator
myGillespie = TauLeap.Sim(cellTypes)
myGillespie.params = myParam
myParam.Hook(myGillespie)
#Set History
myHist = TauLeapParam.Hist(cellTypes)
myGillespie.history = myHist

timeLimit = 1e4

myGillespie.timeLimit = timeLimit
myGillespie.RECORD_TAU_INFO = 0

myGillespie.n_c = 10
myGillespie.stopAtAppear = 0
myGillespie.epsilon = 0.1

#Param Stuff
myParam.n0[0] = 1e4
s = 0.01
for i in range(0,cellTypes):
    myParam.r[i] = math.pow(1.0 + s, i)
    myParam.u[i] = 1e-4
myParam.USE_D = True
myParam.d = 3

#Deter Stuff
deterSim = DeterministicTauLeap.Sim(cellTypes)
deterSim.params = myParam
myHist2 = TauLeapParam.Hist(cellTypes)
deterSim.history = myHist2

deterSim.timeStep = 0.01
deterSim.timeLimit = timeLimit

SDP = 4
PointCount = 8

dataX = []
dataY = []
dataY_deter = []
dataY_deter = []

minN = 1e1
maxN = 1e8

for p in range(0,PointCount):
    startTime = time.clock()
    print("Current Data Point = {0}/{1} ({2}%)".format(p + 1, PointCount, 100.0 * float(p)/(PointCount-1)))
    
    myParam.n0[0] = int(SimUtil.SweepParameterLog(p,PointCount, minN, maxN))
    myParam.u = [1e-2 / myParam.n0[0]] * cellTypes
    print(myParam.u)
    print("N", myParam.n0[0])
    res = myGillespie.SimulateBatch(SDP)
    dataX.append(myParam.n0[0])
    dataY.append(res.avgFixTime)
    
    deterSim.Integrate()  
    dataY_deter.append(deterSim.curTime)

plt.plot(dataX, dataY, 'o')
plt.plot(dataX, dataY_deter)
plt.xscale("log")