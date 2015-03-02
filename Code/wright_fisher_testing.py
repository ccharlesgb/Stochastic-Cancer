# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 18:53:15 2015

@author: Jonny
"""

import wright_fisher
import matplotlib.pyplot as plt
import math

cellTypes = 21
population = 1e9

myWF = wright_fisher.wright_fisher(cellTypes)
myHist = wright_fisher.wf_hist(cellTypes)

myWF.u = 1e-9
#myWF.u = 1e-2 / population

myWF.iN[0] = population

myWF.history = myHist

myWF.popSize = population

myWF.stepLimit = 10000

myWF.useApproxTheta = 0

myWF.d = 200

s = 0.01
adv = 0
for i in range(0,cellTypes):
    if i % 2 == 0:
        myWF.r[i] = math.pow(1.0 + s, adv)
        adv += 1
    else:
        myWF.r[i] = math.pow(1.0 + s, adv)
        adv += 1


myWF.Simulate()

plt.figure()
plt.subplot(211)
for t in range(0, myWF.curStep, myWF.curStep / 8):
    dataX = []
    dataY = []
    for i in range(0, cellTypes):
        dataX.append(i)
        dataY.append(myHist.histArray[i][t])
    plt.plot(dataX,dataY, 'o-')
    
plt.yscale("log")
plt.xlabel("Number of Mutations")
plt.ylabel("Cell count")
plt.show()

plt.subplot(212)
for i in range(0, cellTypes):
   plt.plot(myHist.stepHist, myHist.histArray[i])
   plt.yscale("log")

plt.yscale("log")
plt.show()

plt.figure()
plt.subplot(211)
for i in range(0, cellTypes):
   plt.plot(myHist.stepHist, myHist.thetajHist[i])

plt.subplot(212)
plt.plot(myHist.avgJHist)


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
