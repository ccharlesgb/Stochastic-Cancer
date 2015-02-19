# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 13:45:50 2015

@author: Connor
"""

import matplotlib.pyplot as plt
import GillespieTauLeap
import CarcinogenNParam_OPTON
import math

TYPE_COUNT = 10

myParam = CarcinogenNParam_OPTON.CarcinogenNParam(TYPE_COUNT)

#Create the simulator
myGillespie = GillespieTauLeap.Gillespie()
myGillespie.Hook(myParam) #VERY IMPORTANT

#Set History
myHist = CarcinogenNParam_OPTON.CarcinogenNHist(TYPE_COUNT)
myGillespie.history = myHist

minr1 = 0.3
maxr1 = 3.0

DPC = 10
SDP = 500

dataX = []
dataY = []

myGillespie.timeLimit = 10000

myGillespie.RECORD_TAU_INFO = 1

myParam.n0[0] = 1e6

for i in range(0,TYPE_COUNT):
    myParam.r[i] = math.pow(1.0 + 0.01, i)
    myParam.u[i] = 1.0 / myParam.n0[0]
    
myGillespie.epsilon = 0.1

SDP = 10


res = myGillespie.SimulateBatch(SDP)

plt.figure()

print("STEPS", myGillespie.simSteps)
for t in range(0, myGillespie.simSteps, myGillespie.simSteps / 8):
    dataX = []
    dataY = []
    for i in range(0, TYPE_COUNT):
        dataX.append(i)
        dataY.append(myHist.histArray[i][t])
        plt.plot(dataX,dataY, 'o-')
        plt.yscale("log")

plt.show()