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
import NumericalTau_TEST

def GetTau(param):
    s = param.r[1] - param.r[0]
    top = math.pow(math.log(s / (param.u[0] * param.d)),2.0)
    bottom = 2.0 * s * math.log(param.popSize)
    return top/bottom

def GetTau2(param):
    s = param.r[1] - param.r[0]
    logs = math.log(s / (param.u[0] * param.d)*math.sqrt(2.0*math.log(param.popSize)))
    top = math.pow(logs,2.0)
    bottom = 2.0 * s * math.log(param.popSize)
    return top/bottom

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
    
minU = 1e-8
maxU = 1e-5

SDP = 1
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
    
    myParam.u[0] = SimUtil.SweepParameterLog(p,PointCount, minU, maxU)
    print("U", myParam.u[0])
    res = myWF.SimulateBatch(SDP)
    dataX.append(myParam.u[0])
    dataY.append(res.avgFixTime)
    
    theory = 0.0
    theory2 = 0.0
    
    theory = myParam.AnalyticalWaitingTime()
    for i in range(0, cellTypes + 1):
        theory2 += NumericalTau_TEST.SolveTauIntegral(i, myParam)
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
plt.xlabel("U")
plt.ylabel("t_{0}".format(cellTypes - 1))
plt.xlim(minU, maxU)

plt.subplot(212)
plt.plot(dataX, dataY_anal_err)
plt.plot(dataX, dataY_anal2_err)
plt.xscale("log")
plt.xlabel("U")
plt.ylabel("Error")
plt.xlim(minU, maxU)

plt.show()

plt.figure()
for i in range(0, cellTypes):
   plt.plot(myHist.stepHist[0::1], myHist.histArray[i][0::1])
   plt.yscale("log")

file_name = "WrightFisherSweepU_SDP_{0}_DPC_{1}_CT_{2}_S_{3}_U_{4}".format(SDP,PointCount,cellTypes,s,myParam.u[0])

data = dict()
data["U"] = dataX
data["Ut_20"] = dataY
data["Ut_20_anal"] = dataY_anal
data["Ut_20_anal2"] = dataY_anal2

MatTools.SaveDict(file_name, data)