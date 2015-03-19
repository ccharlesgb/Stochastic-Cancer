# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 16:03:47 2015

@author: Connor
"""

import wright_fisher
import matplotlib.pyplot as plt
import math
import SimUtil
import NumericalTau_TEST
import TauSolver
import numpy as np
import gaussian_fitter

cellTypes = 20
population = 1e9

myWF = wright_fisher.wright_fisher()
myHist = wright_fisher.wf_hist(cellTypes)
myParam = wright_fisher.wright_fisher_params(cellTypes)

myParam.b = 0.01
myParam.d = 100

myWF.stopAtAppear = 1

myParam.iN[0] = population
myWF.history = myHist
myWF.stepLimit = 1000000
myWF.useApproxTheta = 0
myWF.params = myParam


myParam.u[0] = 1e-7

SPD = 1
DPC = 5

s = 0.05
for i in range(0,cellTypes):
    myParam.r[i] = math.pow(1.0 + s, i)
    #myParam.r[i] = 1.0 + s*( (i/(math.sqrt(cellTypes)) )**2 ) 

myWF.reset()

myWF.Simulate()

steps = 10
final_time = myWF.curStep
delta_time = final_time/(steps+2)
average_j = []
average_sj = []
time = []
'''
#helper function to remove zero values from xdata
def trimZeros(xdata, ydata):
    start = -1 
    end = -1
    for i in range(0, cellTypes):
        if(ydata[i] > 0 and start == -1):
            start = i
        if(ydata[i] == 0 and start != -1 ):
            end = i-1
    
    return xdata[start:end]
 '''

plt.figure()
plt.subplot(2,1,1)
for t in range(delta_time, final_time, delta_time):
    time.append(t)    
    dataX = []
    dataY = []    
    for i in range(0, cellTypes):
        dataX.append((i))
        dataY.append(myHist.histArray[i][t])
    plt.plot(dataX, dataY, "-o")
    average_j.append(myHist.avgJHist[t])    
    plt.xlabel("Mutations")
    plt.ylabel("Frequency")   



plt.subplot(2,1,2)
for t in range(delta_time, final_time, delta_time):
    dataX = []
    dataY = []    
    for i in range(0, cellTypes):
        dataX.append(myParam.r[i]) #fitness
        dataY.append(myHist.histArray[i][t])
    plt.plot(dataX, dataY, "-o")
    average_sj.append(myHist.avgSJHist[t])    
    plt.xlabel("Fitness")
    plt.ylabel("Frequency")     
    
plt.xlabel("Time")
plt.ylabel("Frequency")


plt.figure()
plt.subplot(2,1,1)
plt.xlabel("Time")
plt.ylabel("Average j")
plt.plot(time, average_j)

plt.subplot(2,1,2)
plt.xlabel("Time")
plt.ylabel("Average S_j")
plt.plot(time, average_sj)

plt.show()