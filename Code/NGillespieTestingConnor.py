# -*- coding: utf-8 -*-
"""
Created on Tue Feb 03 13:13:05 2015

@author: Connor
"""

import matplotlib.pyplot as plt
import GillespieTauLeap_OPTON
import CarcinogenNParam_OPTON
import math

TYPE_COUNT = 20

myParam = CarcinogenNParam_OPTON.CarcinogenNParam(TYPE_COUNT)
myParam.IJ_DIFF_TOLERANCE = 5

#Create the simulator
myGillespie = GillespieTauLeap_OPTON.Gillespie()
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
    

myGillespie.epsilon = 1.0

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