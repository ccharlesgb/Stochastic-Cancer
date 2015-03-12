# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 18:53:15 2015

@author: Jonny
"""

import wright_fisher
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
population = 1e9

myWF = wright_fisher.wright_fisher()
myHist = wright_fisher.wf_hist(cellTypes)
myParam = wright_fisher.wright_fisher_params(cellTypes)

<<<<<<< HEAD
myParam.u = [1e-7] * 10
<<<<<<< Updated upstream
=======
>>>>>>> origin/master
myParam.d = 100

myWF.stopAtAppear = 1
=======
myParam.d = 200
>>>>>>> Stashed changes

myParam.iN[0] = population
myWF.history = myHist
myWF.stepLimit = 1000000
myWF.useApproxTheta = 0
myWF.params = myParam

<<<<<<< HEAD
<<<<<<< Updated upstream
=======
myParam.uNotConst = 1

>>>>>>> origin/master
s = 0.01
d = 1.8
for i in range(0,cellTypes):
    myParam.r[i] = math.pow(1.0 + s, i)
    #myParam.r[i] = 1.0 + (i*i)*s
<<<<<<< HEAD
=======
s = 0.05
>>>>>>> Stashed changes
=======
    myParam.u[i] = 1e-10 * math.pow(d,i)

print(myParam.u)
>>>>>>> origin/master

for i in range(0,cellTypes):
    #myParam.r[i] = math.pow(1.0 + s, i)
    myParam.r[i] = 1.0 + s*(i**(3))    
    #myParam.r[i] = 1.0 + math.log(1.0 + s)
myWF.Simulate()

plt.figure()
plt.subplot(211)
<<<<<<< Updated upstream
for t in range(0, myWF.curStep, max(myWF.curStep / 8, 1)):
=======
for t in range(0, myWF.curStep, myWF.curStep / 16):
>>>>>>> Stashed changes
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
   plt.plot(myHist.stepHist[0::1], myHist.histArray[i][0::1])
   plt.yscale("log")
plt.xlabel("Time")
plt.ylabel("Cell Count")
plt.yscale("log")
plt.show()

plt.figure()
plt.subplot(311)
for i in range(0, cellTypes):
   plt.plot(myHist.stepHist, myHist.thetajHist[i])
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

print("Appearance Time: {0}".format(myWF.curStep))

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
