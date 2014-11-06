# -*- coding: utf-8 -*-
"""
Created on Thu Nov 06 13:48:28 2014

@author: Connor
"""

import TwoSpecies
import math

import matplotlib.pyplot as plt

mySim = TwoSpecies.Gillespie(10)
mySim.timeLimit = 10000
mySim.u1 = 0.1

dataPointCount = 50

minu1 = 0.002
maxu1 = 0.05

sdp = 100

dataX = []
dataFix = []

dataTheory = []

mySim.r0 = 1.0
mySim.r1 = 1.0

for i in range(0,dataPointCount):
    fixTime = 0.0
    mySim.u1 = float(maxu1 - minu1) * float(i) / (dataPointCount - 1.0) + minu1
    print(mySim.u1)
    for sim in range(0, sdp):
        mySim.Simulate()
        if mySim.Fixated():
            fixTime += mySim.curTime
        else:
            print("WARNING DIDNT REACH FIXATION")
    
    dataX.append(mySim.u1)
    dataFix.append(fixTime / sdp)
    
    rho = 1.0 - (mySim.r0 * (1.0 - mySim.u1) / (mySim.r1 + mySim.r0 * mySim.u1))
    rho = rho / (1.0 - math.pow(mySim.r0 * (1.0 - mySim.u1) / (mySim.r1 + mySim.r0 * mySim.u1), mySim.N))
    Rate = mySim.N * mySim.u1 * rho
    theory = mySim.N / Rate
    dataTheory.append(theory)
    
plt.plot(dataX,dataFix)
plt.plot(dataX,dataTheory)