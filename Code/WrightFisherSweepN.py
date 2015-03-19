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

mySolver = TauSolver.Solver(myParam)    

minN = 1e6
maxN = 1e9

SDP = 25
PointCount = 8

dataX = []
dataY = []
dataY_anal1 = []
dataY_anal2 = []
dataY_anal3 = []
dataY_anal1_err = []
dataY_anal2_err = []
dataY_anal3_err= []

myWF.reset()
for p in range(0,PointCount):
    startTime = time.clock()
    print("Current Data Point = {0}/{1} ({2}%)".format(p + 1, PointCount, 100.0 * float(p)/(PointCount-1)))
    
    myParam.iN[0] = int(SimUtil.SweepParameterLog(p,PointCount, minN, maxN))
    print("N", myParam.iN[0])
    res = myWF.SimulateBatch(SDP)
    dataX.append(myParam.iN[0])
    dataY.append(res.avgFixTime)
    
    mySolver.CacheX0()    
    
    theory = 0.0
    theory2 = 0.0
    
    theory1 = mySolver.GetWaitingTimeOriginal(cellTypes)
    theory2 = mySolver.GetWaitingTime(cellTypes)
    theory3 = mySolver.GetWaitingTimeNeglect(cellTypes)
    
    dataY_anal1.append(theory1)
    dataY_anal2.append(theory2)
    dataY_anal3.append(theory3)
    
    dataY_anal1_err.append(theory1 - res.avgFixTime)
    dataY_anal2_err.append(theory2 - res.avgFixTime)
    dataY_anal3_err.append(theory3 - res.avgFixTime)
    
    print("Complete (Took {:.{s}f} seconds)".format(time.clock() - startTime, s=1))    

plt.figure()
plt.subplot(211)
plt.plot(dataX,dataY, 'o')
plt.plot(dataX,dataY_anal1, ':')
plt.plot(dataX,dataY_anal2, '-')
plt.plot(dataX,dataY_anal3, '--')
plt.xscale("log")
plt.xlabel("N")
plt.ylabel("t_{0}".format(cellTypes - 1))
plt.xlim(minN, maxN)

plt.subplot(212)
plt.plot(dataX, dataY_anal1_err,':')
plt.plot(dataX, dataY_anal2_err, '-')
plt.plot(dataX, dataY_anal3_err, '--')
plt.xscale("log")
plt.xlabel("N")
plt.ylabel("Error")
plt.xlim(minN, maxN)

plt.show()

file_name = "WrightFisherSweepN_SDP_{0}_DPC_{1}_CT_{2}_S_{3}_U_{4}".format(SDP,PointCount,cellTypes,s,myParam.u[0])

data = dict()
data["N"] = dataX
data["Nt_20"] = dataY
data["Nt_20_anal1"] = dataY_anal1
data["Nt_20_anal2_transient"] = dataY_anal2
data["Nt_20_anal3_neglect"] = dataY_anal3

MatTools.SaveDict(file_name, data)
