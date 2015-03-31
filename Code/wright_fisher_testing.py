# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 18:53:15 2015

@author: Jonny
"""

import wright_fisher
import TauLeapParam
import matplotlib.pyplot as plt
import math

def GetTauJ(param):
    summation = 0.0
    s = param.r[1] - param.r[0]
    for i in range(0, param.cellTypes):
        top = math.log(s / (param.u[i] * param.d))
        bottom = 2.0 * s * math.log(param.N)
        summation += (top / bottom)
    return summation;


cellTypes = 21
population = 1e6

myWF = wright_fisher.wright_fisher(cellTypes)
myHist = TauLeapParam.Hist(cellTypes)
myParam = TauLeapParam.Params(cellTypes)

myParam.u = [1e-7] * cellTypes

myParam.d = 100

myWF.stopAtAppear = 1


myParam.n0[0] = population
myWF.history = myHist
myWF.stepLimit = 1000000
myWF.useApproxTheta = 0
myWF.params = myParam


myParam.uNotConst = 0

s = 0.01
d = 1.8
for i in range(0,cellTypes):
    myParam.r[i] = math.pow(1.0 + s, i)
    #myParam.r[i] = 1.0 + (i*i)*s


print(myParam.u)

myWF.Simulate()

plt.figure()
plt.subplot(211)
for t in range(0, myWF.curTime, max(myWF.curTime / 8, 1)):
    dataX = []
    dataY = []
    for i in range(0, cellTypes):
        #dataX.append(myParam.r[i]-1.0)
        dataX.append(i)        
        dataY.append(myHist.histArray[i][t])
    plt.plot(dataX,dataY, 'o-')
    
plt.yscale("log")
plt.xlabel("S_j")
plt.ylabel("Cell count")
plt.show()

plt.subplot(212)
for i in range(0, cellTypes):
   plt.plot(myHist.tHist[0::1], myHist.histArray[i][0::1])
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

'''

s_min = 0.01
s_max = 0.50
steps_to_fix = []
SPD = 10
DPC = 4
s = []

for curPoint in range(0,DPC):
    s.append( s_max*(curPoint/(DPC-1) ) )
    myWF.s = s[curPoint]
    
    fix_step_term = 0    
    for dataPoint in range(0,SPD):
        myWF.Simulate()
        fix_step_term += myWF.curStep
    steps_to_fix.append(float(fix_step_term)/SPD)
    
plt.plot(s,steps_to_fix)
plt.xlabel("Selective advantage")
plt.ylabel("Steps to fixation")
plt.xscale("log")
'''
