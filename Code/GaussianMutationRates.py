# -*- coding: utf-8 -*-
"""
Created on Sat May 02 15:52:00 2015

@author: Connor
"""


import wright_fisher
import matplotlib.pyplot as plt
import math
import SimUtil
import TauSolver
import TauLeapParam
import MatTools
import random
import numpy


def LogNormalU(mean,variance, count):
    u = numpy.random.lognormal(mean,variance, count)
    return u
    
def LogNormalU2(mean,variance,count):
    log_mean = math.log10(mean)
    u = []
    for i in range(0,count):
        log_u = random.normalvariate(log_mean, variance)
        u.append(math.pow(10.0, log_u))
    return u

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

SPD = 100
DPC = 8

myHist.SPD = SPD

s = 1e-2
myParam.SetCompoundFitness(s)   
myParam.SetUAll(1e-4)

minVariance = 0.1
maxVariance = 1

meanU = 1e-6

dataVar = []
dataApp = []

dataAppTimes = []
for i in range(0, DPC):
    dataAppTimes.append([])
    
dataVariances = []

#dataApp_neglect = []

mySolver = TauSolver.Solver(myParam)
mySolver.CacheX0()

for dp in range(0, DPC):
    variance = SimUtil.SweepParameter(dp,DPC, minVariance,maxVariance)
    print("Variance is {0}".format(variance))
    
    for i in range(0,SPD):
        u_arr = LogNormalU2(meanU, variance, cellTypes-1)
        for j in range(0, cellTypes-1):
            myParam.u[j] = u_arr[j]
        print("MYPARAM U ", min(myParam.u), max(myParam.u))
        myWF.Simulate()
        dataAppTimes[dp].append(myWF.curTime)
        
    dataVar.append(variance)
    dataVariances.append(numpy.var(dataAppTimes[dp]))

    #dataApp.append(res.avgFixTime)
    #dataApp_neglect.append(mySolver.GetWaitingTimeNeglect(cellTypes - 1))

plt.figure()
plt.plot(dataVar, dataVariances)
#plt.plot(dataVar, dataApp_neglect)