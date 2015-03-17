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

cellTypes = 21
population = 1e9

myWF = wright_fisher.wright_fisher()
myHist = wright_fisher.wf_hist(cellTypes)
myParam = wright_fisher.wright_fisher_params(cellTypes)

myParam.d = 100

myWF.stopAtAppear = 1

myParam.iN[0] = population
myWF.history = myHist
myWF.stepLimit = 1000000
myWF.useApproxTheta = 0
myWF.params = myParam

myParam.u[0] = 1e-7

s = 0.01
for i in range(0,cellTypes):
    myParam.r[i] = math.pow(1.0 + s, i)

mySolver = TauSolver.Solver(myParam)
    
minS = 1e-4
maxS = 1e-1

SDP = 10
PointCount = 8

dataX = []
dataY = []
dataY_anal = []
dataY_anal2 = []
dataY_anal_err = []
dataY_anal2_err = []

myWF.reset()
for p in range(0,PointCount):
    startTime = time.clock()
    print("Current Data Point = {0}/{1} ({2}%)".format(p + 1, PointCount, 100.0 * float(p)/(PointCount-1)))
    
    s = SimUtil.SweepParameterLog(p,PointCount, minS, maxS)
    for i in range(0,cellTypes):
        myParam.r[i] = math.pow(1.0 + s, i)
    print("S", s)
    res = myWF.SimulateBatch(SDP)
    dataX.append(s)
    dataY.append(res.avgFixTime)
    
    mySolver.CacheX0()
    
    theory = 0.0
    theory2 = 0.0
    
    theory = mySolver.GetWaitingTime(cellTypes)
    theory2 = mySolver.GetWaitingTimeNeglect(cellTypes)
    dataY_anal.append(theory)
    dataY_anal2.append(theory2)
    
    dataY_anal_err.append(theory - res.avgFixTime)
    dataY_anal2_err.append(theory2 - res.avgFixTime)
    
    print("Complete (Took {:.{s}f} seconds)".format(time.clock() - startTime, s=1))    

plt.figure()
plt.subplot(211)
plt.plot(dataX,dataY, 'o')
plt.plot(dataX,dataY_anal, ':')
plt.plot(dataX,dataY_anal2, '--')
plt.xscale("log")
plt.xlabel("S")
plt.ylabel("t_{0}".format(cellTypes - 1))
plt.xlim(minS, maxS)

plt.subplot(212)
plt.plot(dataX, dataY_anal_err)
plt.plot(dataX, dataY_anal2_err)
plt.xscale("log")
plt.xlabel("S")
plt.ylabel("Error")
plt.xlim(minS, maxS)

plt.show()

file_name = "WrightFisherSweepS_SDP_{0}_DPC_{1}_CT_{2}_S_{3}_U_{4}".format(SDP,PointCount,cellTypes,s,myParam.u[0])

data = dict()
data["S"] = dataX
data["St_20"] = dataY
data["St_20_anal"] = dataY_anal
data["St_20_anal_new"] = dataY_anal2

MatTools.SaveDict(file_name, data)