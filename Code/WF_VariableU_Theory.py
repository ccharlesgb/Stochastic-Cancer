# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 18:53:15 2015

@author: Jonny
"""

import wright_fisher
import matplotlib.pyplot as plt
import SimUtil
import math
import TauLeapParam

def GetT_K(param):
    summation = 0.0
    s = param.r[1] - param.r[0]
    for i in range(0, param.typeCount - 1):
        top = math.pow(math.log(s / (param.u[i] * param.d)), 2.0)
        bottom = 2.0 * s * math.log(param.N)
        summation += (top / bottom)
    return summation;


cellTypes = 21
population = 1e9

myWF = wright_fisher.wright_fisher(cellTypes)
myHist = TauLeapParam.Hist(cellTypes)
myParam = TauLeapParam.Params(cellTypes)

myParam.d = 100

myWF.stopAtAppear = 1

myParam.n0[0] = population
#myWF.history = myHist
myWF.stepLimit = 1000000
myWF.useApproxTheta = 0
myWF.params = myParam
myWF.printProgress = 1

myParam.uNotConst = 1

s = 0.01
uPow = 2.0

uPowMin = 1.8
uPowMax = 2.2

SPD = 5
DPC = 5

dataX = []
dataY = []
dataY_Theory = []

s = 0.01
for i in range(0,cellTypes):
    myParam.r[i] = math.pow(1.0 + s, i)

for curPoint in range(0,DPC):
    uPow = SimUtil.SweepParameter(curPoint, DPC, uPowMin, uPowMax)
    for i in range(0,cellTypes-1):
        myParam.u[i] = 1e-10 * math.pow(uPow,i)
    
    print("U: ", myParam.u[0], myParam.u[cellTypes - 2])    

    dataX.append(uPow)

    res = myWF.SimulateBatch(SPD)
    dataY.append(res.avgFixTime)
    dataY_Theory.append(GetT_K(myParam))
    
plt.plot(dataX, dataY, 'o')
plt.plot(dataX, dataY_Theory)
plt.xlabel("uPow")
plt.ylabel("Steps to fixation")

