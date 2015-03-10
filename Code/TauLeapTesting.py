# -*- coding: utf-8 -*-
"""
Created on Mon Mar 09 15:09:08 2015

@author: Connor
"""

import matplotlib.pyplot as plt
import MatTools
import numpy as np
import math
import time

import TauLeap
import TauLeapParam

cellTypes = 21
myParam = TauLeapParam.Params(cellTypes)

#Create the simulator
myGillespie = TauLeap.Sim(cellTypes)
myGillespie.params = myParam
myParam.Hook(myGillespie)
#Set History
myHist = TauLeapParam.Hist(cellTypes)
myGillespie.history = myHist

minr1 = 0.3
maxr1 = 3.0

DPC = 10
SDP = 1

dataX = []
dataY = []

myGillespie.timeLimit = 100000

myGillespie.RECORD_TAU_INFO = 1

myParam.n0[0] = 1e9

s = 0.01
for i in range(0,cellTypes):
    myParam.r[i] = math.pow(1.0 + s, i)
    myParam.u[i] = 1e-7

myParam.USE_D = True
myParam.d = 100

myGillespie.n_c = 10
myGillespie.stopAtAppear = 1
myGillespie.epsilon = 0.2

'''
for i in range(0, DPC):
    myParam.r1 = ((float(i) / (DPC - 1)) * (maxr1 - minr1)) + minr1
    dataX.append(myParam.r1)
    print(myParam.r1)
    
    res = myGillespie.SimulateBatch(SDP)

    print("AVERAGE BAD FRAMES: ", float(res.avgBadFrames), float(res.avgFrames))
    
    dataY.append(res.avgFixProb)

plt.plot(dataX, dataY, linewidth=1.0, label="X2(t)")
'''

avgBadFrame = 0
TAU_HIST_TOTAL = []

startTime = time.clock()

for i in range(0,SDP):
    myGillespie.Simulate()
    print("AVERAGE BAD FRAMES: ", float(myGillespie.BAD_FRAME_COUNT))
    avgBadFrame += myGillespie.BAD_FRAME_COUNT
    TAU_HIST_TOTAL.extend(myGillespie.TAU_HIST)

print("Complete (Took {:.{s}f} seconds)".format((time.clock() - startTime)/SDP, s=2))    

plt.figure()
for i in range(0, cellTypes):
   plt.plot(myHist.tHist[0::1], myHist.histArray[i][0::1])
   plt.yscale("log")

if myGillespie.RECORD_TAU_INFO:
    plt.figure()
    weights = np.ones_like(TAU_HIST_TOTAL)/len(TAU_HIST_TOTAL)
    plt.hist(TAU_HIST_TOTAL)
    plt.yscale('log', nonposy='clip')
    
file_name = "tauHistory_epsilon_{0}".format(myGillespie.epsilon)
data = dict()
data["tau"] = TAU_HIST_TOTAL
data["bad_frame"] = avgBadFrame
data["epsilon"] = myGillespie.epsilon

MatTools.SaveDict(file_name,data)
