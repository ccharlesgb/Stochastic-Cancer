# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 16:23:16 2015

@author: Connor
"""

import matplotlib.pyplot as plt
import numpy as np
import math
import time

import MatTools
import DeterministicTauLeap
import TauLeap
import TauLeapParam
import SimUtil
import wright_fisher

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

myWF = wright_fisher.wright_fisher(cellTypes)
myWF.params = myParam

#Param Stuff
myParam.n0[0] = 1e2
s = 0.0
for i in range(0,cellTypes):
    myParam.r[i] = math.pow(1.0 + s, i)

myParam.SetUAll(1e-4)
myParam.USE_D = False
myParam.d = 100

myParam.r[0] = 1.0
myParam.r[1] = 1.0
myParam.r[2] = 1.0

myParam.u[0] = 0.1
myParam.u[1] = 0.1

#Deter Stuff
deterSim = DeterministicTauLeap.Sim(cellTypes)
deterSim.params = myParam
myHist2 = TauLeapParam.Hist(cellTypes)
deterSim.history = myHist2

deterSim.timeStep = 0.01
deterSim.timeLimit = timeLimit
deterSim.stopAtAppear = 0

'''
deterSim.Integrate()
myGillespie.Simulate()

for i in range(0, cellTypes):
    plt.plot(myHist2.tHist, myHist2.histArray[i])
    plt.plot(myHist.tHist, myHist.histArray[i])
'''

SDP = 1
PointCount = 8

dataX = []
dataY = []
dataY_deter = []
dataY_WF = []

minN = 1e1
maxN = 1e2

for p in range(0,PointCount):
    startTime = time.clock()
    print("Current Data Point = {0}/{1} ({2}%)".format(p + 1, PointCount, 100.0 * float(p)/(PointCount-1)))
    
    #myParam.n0[0] = int(SimUtil.SweepParameterLog(p,PointCount, minN, maxN))
    myParam.n0[0] = int(SimUtil.SweepParameter(p,PointCount, minN, maxN))
    #myParam.SetUAll([0.1] * cellTypes
    print(myParam.u)
    print("N", myParam.n0[0])
    res = myGillespie.SimulateBatch(SDP)
    dataX.append(myParam.n0[0])
    dataY.append(res.avgFixTime)
    
    myWF.SimulateBatch(SDP)
    dataY_WF.append(res.avgFixTime)
    
    deterSim.Integrate()  
    dataY_deter.append(deterSim.curTime)

plt.plot(dataX, dataY, 'o')
plt.plot(dataX, dataY_deter)
plt.plot(dataX, dataY_WF)
#plt.xscale("log")

file_name = "TauLeapSweepN_SDP_{0}_DPC_{1}".format(SDP, PointCount)

data = dict()
data["N"] = dataX
data["FixTime"] = dataY

MatTools.SaveDict(file_name, data)

