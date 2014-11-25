# -*- coding: utf-8 -*-
"""
Created on Thu Nov 06 13:48:28 2014

@author: Connor
"""

import TwoSpecies
import math
import anfixtime

import FixationTimeOnePulse2

import matplotlib.pyplot as plt

mySim = TwoSpecies.Gillespie(10)
mySim.timeLimit = 100000
mySim.u1 = 0.1
mySim.ij=0
dataPointCount = 20

minu1 = 0.002
maxu1 = 0.05

sdp = 200 #simsperdatapoint

dataX = []
dataFix = []

dataTheory = []

mySim.r0 = 1.5
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
    
    theory = anfixtime.GetFixTimeJ(mySim, mySim.ij) 
    dataTheory.append(theory)
    
    sumwk = 0.0
    for k in range(1, mySim.N + 1):
        sumwk += FixationTimeOnePulse2.Getwk(mySim, k)
    
    #dataTheory.append(sumwk)

plt.plot(dataX,dataFix, '^')
plt.plot(dataX,dataTheory)