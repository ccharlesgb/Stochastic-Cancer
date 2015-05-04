# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 18:53:15 2015

@author: Jonny
"""

import wright_fisher
import TauLeapParam
import matplotlib.pyplot as plt
import math
import MatTools

def GetTauJ(param):
    summation = 0.0
    s = param.r[1] - param.r[0]
    for i in range(0, param.cellTypes):
        top = math.log(s / (param.u[i] * param.d))
        bottom = 2.0 * s * math.log(param.N)
        summation += (top / bottom)
    return summation;


cellTypes = 21
population = 1e9

myWF = wright_fisher.wright_fisher(cellTypes)
myHist = TauLeapParam.Hist(cellTypes)
myParam = TauLeapParam.Params(cellTypes)

myParam.SetUAll(1e-7)

myParam.d = 100

myWF.stopAtAppear = 1


myParam.n0[0] = population
myWF.history = myHist
myWF.stepLimit = 1000000
myWF.useApproxTheta = 0
myWF.params = myParam


myParam.uNotConst = 0

s = 0.01
myParam.SetCompoundFitness(s)

print(myParam.u)

myWF.Simulate()

save_data = dict()

plt.figure()
plt.subplot(211)
for t in range(0, myWF.curTime, 600):
    dataX = []
    dataY = []
    for i in range(0, cellTypes):
        #dataX.append(myParam.r[i]-1.0)
        dataX.append(i)        
        dataY.append(myHist.histArray[i][t])
    save_data["t_{0}_j".format(t)] = dataX
    save_data["t_{0}_x_j".format(t)] = dataY
    plt.plot(dataX,dataY, 'o-')
    
plt.yscale("log")
plt.xlabel("S_j")
plt.ylabel("Cell count")
plt.show()

plt.subplot(212)
for i in range(0, cellTypes):
   plt.plot(myHist.tHist, myHist.histArray[i])
   plt.yscale("log")
plt.xlabel("Time")
plt.ylabel("Cell Count")
plt.yscale("log")
plt.show()

plt.figure()
plt.subplot(311)
for i in range(0, cellTypes):
   plt.plot(myHist.tHist, myHist.thetajHist[i])
plt.xlabel("Time")
plt.ylabel("Theta J")

plt.subplot(312)
plt.plot(myHist.avgJHist)
plt.xlabel("Time")
plt.ylabel("Average S_j")

gradJ = []
for i in range(500, len(myHist.avgJHist)):
    gradJ.append(myHist.avgJHist[i] - myHist.avgJHist[i-500])
plt.subplot(313)
plt.plot(gradJ)

print("Appearance Time: {0}".format(myWF.curTime))

save_data["hist_t"] = myHist.tHist
for i in range(0, cellTypes):
    save_data["hist_n{0}".format(i)] = myHist.histArray[i]

MatTools.SaveDict2(save_data, PARAMS = myParam.GetFileString())
