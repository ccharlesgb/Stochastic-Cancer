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

cellTypes = 21
population = 1e9

myWF = wright_fisher.wright_fisher(cellTypes)
myHist = wright_fisher.wf_hist(cellTypes)

myWF.u = 1e-7

myWF.popSize = population
myWF.history = myHist

myWF.stepLimit = 10000

s = 0.01
myWF.SetCompoundFitness(s)


minN = 1e6
maxN = 1e9

SDP = 10
PointCount = 4

dataX = []
dataY = []
dataY_anal = []

for p in range(0,PointCount):
    startTime = time.clock()
    print("Current Data Point = {0}/{1} ({2}%)".format(p + 1, PointCount, 100.0 * float(p)/(PointCount-1)))
    
    myWF.iN[0] = round(SimUtil.SweepParameterLog(p,PointCount, minN, maxN))
    print("IN0", myWF.iN[0])
    res = myWF.SimulateBatch(SDP)
    dataX.append(myWF.iN[0])
    dataY.append(res.avgFixTime)
    
    dataY_anal.append(myWF.AnalyticalWaitingTime())
    
    print("Complete (Took {:.{s}f} seconds)".format(time.clock() - startTime, s=1))    
    
plt.plot(dataX,dataY, 'o')
plt.plot(dataX,dataY_anal)
plt.xscale("log")
plt.xlabel("N")
plt.ylabel("t_{0}".format(cellTypes - 1))
plt.xlim(minN, maxN)
plt.show()

file_name = "WrightFisherSweepN_SDP_{0}_DPC_{1}_CT_{2}_S_{3}_U_{4}".format(SDP,PointCount,cellTypes,s,myWF.u)

data = dict()
data["N"] = dataX
data["Nt_20"] = dataY
data["Nt_20_anal"] = dataY_anal

MatTools.SaveDict(file_name, data)