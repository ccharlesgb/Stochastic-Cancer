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
myGillespie.history = myHist

minr1 = 0.3
maxr1 = 3.0

DPC = 10
SDP = 500

dataX = []
dataY = []

myGillespie.timeLimit = 100

myParam.addRate = 0.0

myParam.n0[0] = 10
myParam.n0[1] = 0
myParam.n0[2] = 0

myParam.c0 = 0.001
myParam.c1 = 0.001
myParam.c2 = 0.001

myParam.u1 = 0.1
myParam.u2 = 0.1

myGillespie.epsilon = 1.0

myGillespie.Simulate()


for i in range(0, DPC):
    myParam.r1 = ((float(i) / (DPC - 1)) * (maxr1 - minr1)) + minr1
    dataX.append(myParam.r1)
    print(myParam.r1)
    fixCount = 0
    
    avgBadFrame = 0
    avgTotFrame = 0
    
    for i in range(0,SDP):
        myGillespie.Simulate()
        avgBadFrame += myGillespie.BAD_FRAME_COUNT
        avgTotFrame += myGillespie.simSteps
        if myParam.n[2] == myParam.N:
            fixCount += 1

    print("AVERAGE BAD FRAMES: ", float(avgBadFrame) / SDP, float(avgTotFrame) / SDP)    
    
    fixProb = float(fixCount) / SDP 
    dataY.append(fixProb)

plt.plot(dataX, dataY, linewidth=1.0, label="X2(t)")

plt.figure()
plt.hist(myGillespie.TAU_HIST)

'''

plt.plot(myHist.tHist, myHist.n0Hist, 'o')
plt.plot(myHist.tHist, myHist.n1Hist)
plt.plot(myHist.tHist, myHist.n2Hist)
'''