# -*- coding: utf-8 -*-
"""
Created on Tue Mar 03 14:11:02 2015

@author: Connor
"""

import wright_fisher
import MatTools
import SimUtil
import math
import matplotlib.pyplot as plt

CELL_TYPES = 21
population = 1e9

myParam = wright_fisher.wright_fisher_params(CELL_TYPES)

myParam.u = [1e-7] * 10
myParam.d = 100

myParam.iN[0] = population

xj = []

xjHist = []
avgJHist = []

s = 0.01
for i in range(0,CELL_TYPES):
    myParam.r[i] = math.pow(1.0 + s, i)
    xj.append(0.0)
    xjHist.append([])
    
xj[0] = 1.0

myParam.Reset()

deltaT = 1.0
timeLimit = 6000

curTime = 0.0
for step in range(0, int(timeLimit / deltaT)):
    avgJ = 0.0
    for i in range(0, CELL_TYPES):
        deriv = myParam.GetXJDot(i, xj)
        xj[i] = xj[i] + deriv * deltaT
        curTime = step * deltaT 
        xjHist[i].append(round(xj[i] * population))
        avgJ += i * xj[i]
    avgJHist.append(avgJ)
    if xj[CELL_TYPES - 1] >= 1.0 / population:
        break

plt.figure()
plt.subplot(211)
for t in range(0, step, step / 8):
    dataX = []
    dataY = []
    for i in range(0, CELL_TYPES):
        dataX.append(i)
        dataY.append(xjHist[i][t])
    plt.plot(dataX,dataY, 'o-')
    
plt.yscale("log")
plt.xlabel("Number of Mutations")
plt.ylabel("Cell count")

plt.subplot(212)
for i in range(0, CELL_TYPES):
   plt.plot(xjHist[i])
   plt.yscale("log")

plt.yscale("log")
plt.show()

plt.figure()
plt.plot(avgJHist)
