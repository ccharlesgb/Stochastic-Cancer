# -*- coding: utf-8 -*-
"""
Created on Tue Feb 03 13:13:05 2015

@author: Connor
"""

import matplotlib.pyplot as plt
import GillespieTauLeap
import CarcinogenNParam
import math

TYPE_COUNT = 10

myParam = CarcinogenNParam.CarcinogenNParam(TYPE_COUNT)

#Create the simulator
myGillespie = GillespieTauLeap.Gillespie()
myGillespie.Hook(myParam) #VERY IMPORTANT

#Set History
myHist = CarcinogenNParam.CarcinogenNHist(TYPE_COUNT)
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

myGillespie.Simulate()

plt.figure()
for i in range(0, TYPE_COUNT):
   plt.plot(myHist.tHist, myHist.histArray[i])
   plt.yscale("log")

plt.show()


plt.figure()
plt.subplot(211)
plt.hist(myGillespie.TAU_HIST, log=1)
plt.xlabel("Tau")
plt.ylabel("Frequency")

plt.subplot(212)
plt.plot(myGillespie.TAU_HIST)
plt.xlabel("Step")
plt.ylabel("Tau")
plt.show()