# -*- coding: utf-8 -*-
"""
Created on Sun Apr 12 17:34:59 2015

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
import random

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
myParam.uNotConst = True
s = 0.01
for i in range(0,cellTypes):
    myParam.r[i] = math.pow(1.0 + s, i)

mySolver = TauSolver.Solver(myParam)    

SDP = 50
PointCount = 4

dataX = []
dataY = []
dataY_rand = []

min_var = 0.5
max_var = 3

mean_u = -6

for p in range(0,PointCount):
    startTime = time.clock()
    print("Current Data Point = {0}/{1} ({2}%)".format(p + 1, PointCount, 100.0 * float(p)/(PointCount-1)))
    

    var_u = SimUtil.SweepParameterLog(p, PointCount, min_var, max_var)
    dataX.append(var_u)
    
    print("Start random")
    fixTimeTotal = 0.0
    for i in range(0, SDP):
        for i in range(1,cellTypes):
            ran_u = random.gauss(mean_u, var_u)
            #ran_u = max(ran_u, 1e-10)
            myParam.SetU(i, math.pow(10.0,ran_u))
        print(min(myParam.u), max(myParam.u))
        myWF.Simulate()
        fixTimeTotal += myWF.curTime
    dataY_rand.append(float(fixTimeTotal) / SDP)
    print("Start mean")
    myParam.SetUAll(math.pow(10.0,mean_u))
    res = myWF.SimulateBatch(SDP)
    dataY.append(res.avgFixTime)
    
    print("Complete (Took {:.{s}f} seconds)".format(time.clock() - startTime, s=1))    

plt.figure()
plt.plot(dataX,dataY, 'o-')
plt.plot(dataX,dataY_rand, 'x-')
plt.xlabel("Var")
plt.ylabel("t_{0}".format(cellTypes - 1))

plt.show()

data = dict()


MatTools.SaveDict2(data,SDP = SDP, DPC = PointCount, PARAM = myParam.GetFileString())
