# -*- coding: utf-8 -*-
"""
Created on Thu Nov 06 13:48:28 2014

@author: Connor
"""

import TwoSpecies
import math
import anfixtime
import MatTools


#import FixationTimeOnePulse2

import matplotlib.pyplot as plt

mySim = TwoSpecies.Gillespie(10)
mySim.timeLimit = 1000000
mySim.u1 = 0.1
mySim.ij=0
dataPointCount = 10

mySim.r0 = 1.0
mySim.r1 = 1.0

mySim.u1 = 0.1
mySim.u2 = 0.1

minu1 = 10
maxu1 = 100

sdp = 10000 #simsperdatapoint

dataX = []
dataFix = []

dataTheory = []

mySim.r0 = 1.0
mySim.r1 = 1.0

for i in range(0,dataPointCount):
    fixTime = 0.0
    #mySim.r1 = float(maxu1 - minu1) * float(i) / (dataPointCount - 1.0) + minu1
    mySim.N = int(float(maxu1 - minu1) * float(i) / (dataPointCount - 1.0) + minu1)
    print(mySim.N)
    for sim in range(0, sdp):
        mySim.Simulate()
        if mySim.Fixated():
            fixTime += mySim.curTime
        else:            
            print("WARNING DIDNT REACH FIXATION")
    
    dataX.append(mySim.N)
    dataFix.append(fixTime / sdp)
    
    theory = anfixtime.GetFixTimeJ(mySim, mySim.ij) 
    dataTheory.append(theory)
    
    sumwk = 0.0
    #for k in range(1, mySim.N + 1):
    #    sumwk += FixationTimeOnePulse2.Getwk(mySim, k)
    
    #dataTheory.append(sumwk)

plt.plot(dataX,dataFix, '^')
plt.plot(dataX,dataTheory)

file_name = "TwoSpeciesFix_Anal_sweepN_N_{0}_SDP_{1}".format(mySim.N, sdp)

data2 = dict()
data2["num_fixtime"] = dataTheory

MatTools.SaveXYData(file_name, dataX, dataFix, xLabel = "N", yLabel = "fixtime", sim = mySim, otherDict = data2)