# -*- coding: utf-8 -*-
"""
Created on Tue Feb 03 13:13:05 2015

@author: Connor
"""

import matplotlib.pyplot as plt
import GillespieTauLeap
import CarcinogenParam

myParam = CarcinogenParam.CarcinogenParam()

#Create the simulator
myGillespie = GillespieTauLeap.Gillespie()
myGillespie.Hook(myParam) #VERY IMPORTANT

#Set History
myHist = CarcinogenParam.CarcinogenHist()
#myGillespie.history = myHist

minr1 = 0.3
maxr1 = 3.0

DPC = 10
SDP = 500

dataX = []
dataY = []

myGillespie.timeLimit = 100

myGillespie.RECORD_TAU_INFO = 1

myParam.addRate = 0.0

myParam.n0[0] = 100
myParam.n0[1] = 0
myParam.n0[2] = 0

myParam.c0 = 0.001
myParam.c1 = 0.001
myParam.c2 = 0.001

myParam.u1 = 0.1
myParam.u2 = 0.1

myGillespie.epsilon = 0.1

myGillespie.Simulate()

for i in range(0, DPC):
    myParam.r1 = ((float(i) / (DPC - 1)) * (maxr1 - minr1)) + minr1
    dataX.append(myParam.r1)
    print(myParam.r1)
    
    res = myGillespie.SimulateBatch(SDP)

    print("AVERAGE BAD FRAMES: ", float(res.avgBadFrames), float(res.avgFrames))
    
    dataY.append(res.avgFixProb)

plt.plot(dataX, dataY, linewidth=1.0, label="X2(t)")

'''
plt.plot(myHist.tHist, myHist.n0Hist)
plt.plot(myHist.tHist, myHist.n1Hist)
plt.plot(myHist.tHist, myHist.n2Hist)
'''
if myGillespie.RECORD_TAU_INFO:
    plt.figure()
    plt.hist(myGillespie.TAU_HIST)
