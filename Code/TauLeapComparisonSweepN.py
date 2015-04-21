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
import anfixtime

cellTypes = 2
myParam = TauLeapParam.Params(cellTypes)

#Create the simulator
myGillespie = TauLeap.Sim(cellTypes)
myGillespie.params = myParam
myParam.Hook(myGillespie)

timeLimit = 1e3

myGillespie.timeLimit = timeLimit
myGillespie.RECORD_TAU_INFO = 0

myGillespie.n_c = 10
myGillespie.stopAtAppear = 0
myGillespie.epsilon = 0.05

myWF = wright_fisher.wright_fisher(cellTypes)
myWF.params = myParam
myWF.stopAtAppear = 0

#Set History
myHist = TauLeapParam.Hist(cellTypes)
myWF.history = myHist

myParam.SetUAll(1e-4)
myParam.USE_D = False
myParam.d = 100

myParam.r[0] = 1.0
myParam.r[1] = 1.0
#myParam.r[2] = 1.0

myParam.u[0] = 0.1
#myParam.u[1] = 0.1

SDP = 200
PointCount = 4

dataX = []
dataY = []
dataY_WF = []
data_anal = []

minN = 1e1
maxN = 1e3

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

    #myParam.n0[0] = int(float(myParam.n0[0]) / math.sqrt(2.0))
    res = myWF.SimulateBatch(SDP)
    dataY_WF.append(res.avgFixTime)
    
    myGillespie.Reset()
    #data_anal.append(anfixtime.GetFixTimeJ(myGillespie, myParam,0))

plt.figure()
plt.subplot(211)
plt.plot(myHist.tHist, myHist.histArray[0])
plt.plot(myHist.tHist, myHist.histArray[1])

plt.subplot(212)
plt.plot(dataX, dataY, 'o')
plt.plot(dataX, dataY_WF)
#plt.plot(dataX, data_anal)
#plt.xscale("log")

file_name = "TauLeapSweepN_SDP_{0}_DPC_{1}".format(SDP, PointCount)

data = dict()
data["N"] = dataX
data["FixTime"] = dataY

MatTools.SaveDict(file_name, data)

