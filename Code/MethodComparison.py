# -*- coding: utf-8 -*-
"""
Created on Mon Mar 02 15:35:16 2015

@author: Connor
"""

import matplotlib.pyplot as plt
import GillespieTauLeap_OPTON
import CarcinogenNParam_OPTON
import math
import MatTools

TYPE_COUNT = 11

myParam = CarcinogenNParam_OPTON.CarcinogenNParam(TYPE_COUNT)
myParam.IJ_DIFF_TOLERANCE = 12

#Create the simulator
myGillespie = GillespieTauLeap_OPTON.Gillespie()
myGillespie.Hook(myParam) #VERY IMPORTANT

#Set History
myHist = CarcinogenNParam_OPTON.CarcinogenNHist(TYPE_COUNT)
myGillespie.history = myHist

myGillespie.timeLimit = 10000
myGillespie.RECORD_TAU_INFO = 0
myGillespie.printProgress = 1
myGillespie.stopAtAppear = 1
myParam.n0[0] = 1e9

s = 0.01
for i in range(0,TYPE_COUNT):
    myParam.r[i] = math.pow(1.0 + s, i)
    myParam.u[i] = 1e-7
    
myGillespie.epsilon = 0.1

myGillespie.Simulate()

plt.figure()
for i in range(0, TYPE_COUNT):
   plt.plot(myHist.tHist, myHist.histArray[i])
   plt.yscale("log")

plt.show()

print("Appearance Time: {0}".format(myGillespie.curTime))
