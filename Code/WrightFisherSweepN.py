# -*- coding: utf-8 -*-
"""
Created on Thu Feb 26 13:32:23 2015
@author: Connor
"""

import wright_fisher
import matplotlib.pyplot as plt
import math
import SimUtil
import time
import MatTools
import TauSolver
import TauLeapParam

cellTypes = 21
population = 1e9

myWF = wright_fisher.wright_fisher(cellTypes)
myHist = TauLeapParam.Hist(cellTypes)
myParam = TauLeapParam.Params(cellTypes)

myParam.d = 100

myWF.stopAtAppear = 1

myParam.n0[0] = population
myWF.history = myHist
myWF.timeLimit = 1000000
myWF.useApproxTheta = 0
myWF.params = myParam

myParam.SetUAll(1e-7)

s = 0.01
for i in range(0,cellTypes):
    myParam.r[i] = math.pow(1.0 + s, i)

mySolver = TauSolver.Solver(myParam)    

minN = 1e6
maxN = 1e9

SDP = 25
PointCount = 4

dataX = []
dataY = []
dataY_anal1 = []
dataY_anal2 = []
dataY_anal3 = []
dataY_anal4 = []
dataY_anal1_err = []
dataY_anal2_err = []
dataY_anal3_err = []
dataY_anal4_err = []


for p in range(0,PointCount):
    startTime = time.clock()
    print("Current Data Point = {0}/{1} ({2}%)".format(p + 1, PointCount, 100.0 * float(p)/(PointCount-1)))
    
    myParam.n0[0] = int(SimUtil.SweepParameterLog(p,PointCount, minN, maxN))
    print("N", myParam.n0[0])
    res = myWF.SimulateBatch(SDP)
    dataX.append(myParam.n0[0])
    dataY.append(res.avgFixTime)
    
    mySolver.CacheX0()    
    
    theory = 0.0
    theory2 = 0.0
    theory1 = mySolver.GetWaitingTimeOriginal(cellTypes - 1)
    theory2 = mySolver.GetWaitingTimeRecursive(cellTypes - 1)
    theory3 = mySolver.GetWaitingTimeNeglect(cellTypes - 1)
    theory4 = mySolver.GetWaitingTimeCorrect(cellTypes - 1)
    
    dataY_anal1.append(theory1)
    dataY_anal2.append(theory2)
    dataY_anal3.append(theory3)
    dataY_anal4.append(theory4)
    
    dataY_anal1_err.append(theory1 - res.avgFixTime)
    dataY_anal2_err.append(theory2 - res.avgFixTime)
    dataY_anal3_err.append(theory3 - res.avgFixTime)
    dataY_anal4_err.append(theory4 - res.avgFixTime)
    
    print("Complete (Took {:.{s}f} seconds)".format(time.clock() - startTime, s=1))    

plt.figure()
plt.subplot(211)
plt.plot(dataX,dataY, 'o')
plt.plot(dataX,dataY_anal1, ':')
plt.plot(dataX,dataY_anal2, '-')
plt.plot(dataX,dataY_anal3, '--')
plt.plot(dataX,dataY_anal4, '-.')
plt.xscale("log")
plt.xlabel("N")
plt.ylabel("t_{0}".format(cellTypes - 1))
plt.xlim(minN, maxN)

plt.subplot(212)
plt.plot(dataX, dataY_anal1_err,':')
plt.plot(dataX, dataY_anal2_err, '-')
plt.plot(dataX, dataY_anal3_err, '--')
plt.plot(dataX, dataY_anal4_err, '-.')
plt.xscale("log")
plt.xlabel("N")
plt.ylabel("Error")
plt.xlim(minN, maxN)

plt.show()

data = dict()
data["N"] = dataX
data["Nt_20"] = dataY
data["Nt_20_anal1"] = dataY_anal1
data["Nt_20_anal2_recursive"] = dataY_anal2
data["Nt_20_anal3_neglect"] = dataY_anal3
data["Nt_20_anal4_correct"] = dataY_anal4

MatTools.SaveDict2(data,SDP = SDP, DPC = PointCount, PARAM = myParam.GetFileString())
