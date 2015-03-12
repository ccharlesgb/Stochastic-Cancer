# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 16:03:47 2015

@author: Connor
"""

import wright_fisher
import matplotlib.pyplot as plt
import math
import SimUtil

def GetTau2(param):
    s = param.r[1] - param.r[0]
    logs = math.log(s / (param.u[0] * param.d)*math.sqrt(2.0*math.log(param.popSize)))
    top = math.pow(logs,2.0)
    bottom = 2.0 * s * math.log(param.popSize)
    return top/bottom

def GetTau(param):
    s = param.r[1] - param.r[0]
    top = math.pow(math.log(s / (param.u[0] * param.d)),2.0)
    bottom = 2.0 * s * math.log(param.popSize)
    return top/bottom


cellTypes = 11
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

uPowMin = 1.8
uPowMax = 2.2

SPD = 50
DPC = 5

s = 0.01
for i in range(0,cellTypes):
    myParam.r[i] = math.pow(1.0 + s, i)


avgAppear = [0.0] * cellTypes
theoreticalAppear = [0.0] * cellTypes
theoreticalAppear2 = [0.0] * cellTypes

for curPoint in range(0,SPD):
    myWF.reset()
    myWF.Simulate()       
   
    #Loop through each step
    for i in range(0, cellTypes):
        found = False
        for t in range(len(myHist.histArray[i])):
            if found:
                continue
            count = myHist.histArray[i][t]
            if count > 0:
                avgAppear[i] += t
                found = True
                
for i in range(0, cellTypes):
    theoreticalAppear[i] = i * GetTau(myParam)
    avgAppear[i] /= SPD
    print("Avg Appear {0} = {1}".format(i,avgAppear[i]))
    if (i > 0):
        print("TAU {0} = {1}".format(i, avgAppear[i]-avgAppear[i-1]))
for i in range(2, cellTypes):
    theoreticalAppear2[i] = (i-2) * GetTau2(myParam)

plt.figure()
plt.subplot(211)
plt.plot(avgAppear, 'o')
plt.plot(theoreticalAppear)
plt.plot(theoreticalAppear2, '--')
plt.xlabel("i")
plt.ylabel("Appearance Time")

plt.subplot(212)
for i in range(0, cellTypes):
   plt.plot(myHist.stepHist[0::1], myHist.histArray[i][0::1])
   plt.yscale("log")

plt.yscale("log")
plt.show()

