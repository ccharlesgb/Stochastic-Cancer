# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 14:52:31 2015

@author: Connor
"""

import math
import random
import matplotlib.pyplot as plt
import numpy as np
import MatTools


def GetAvgWaitTime(TYPE_COUNT, SDP, u = 1e-7, s = 1e-2):
    avgWaitTime = 0.0
    for i in range(0, SDP):
        X = []
        pVec = []
    
        CELL_COUNT = int(1e9)
        total = 0.0
        for i in range(0, TYPE_COUNT):
            pVec.append(0.0)
            total += pVec[i]
            X.append(0)
        
        STEPS = 30000
        
        X[0] = CELL_COUNT
        
        r = []
        for i in range(0,TYPE_COUNT):
            r.append(math.pow(1.0 + 0.01, i))
        
        tHist = []
        histArray = dict()
        for i in range(0, TYPE_COUNT):
            histArray[i] = []
        
        step = 0
        fixed = 0
        while (step < STEPS and fixed != 1):
            total = 0.0
            for i in range(0, TYPE_COUNT):
                avgFit = 0.0
                for l in range(0, TYPE_COUNT):
                    avgFit += r[l]*X[l]
                #print(avgFit)
                theta_i = (r[i] * X[i])/avgFit
                if i > 0:
                    theta_i += float(u * (TYPE_COUNT - i + 1) * r[i-1]*X[i-1] )/ avgFit
                    
                pVec[i] = theta_i
                total += pVec[i]
                
            for i in range(0, TYPE_COUNT):
                pVec[i] /= total
                #if pVec[i] == 1.0:
                    #print("FIXATION")
                    #fixed = 1
        
            X = np.random.multinomial(CELL_COUNT, pVec)
            
            if X[TYPE_COUNT - 1] >= 1:
                fixed = 1
                avgWaitTime += step
            
            tHist.append(step)
            for i in range(0, TYPE_COUNT):
                histArray[i].append(X[i])
                
            step += 1
    return float(avgWaitTime) / SDP
    
dataX = []
dataY = []

SDP = 1

minU = 1e-2
maxU = 1e-7

DPC = 2

for i in range(0,DPC):
    print("DATA POINT" , i)
    U = (maxU - minU) * float(i) / (DPC - 1) + minU
    dataX.append(U)
    dataY.append(GetAvgWaitTime(40, SDP, U))
    
plt.plot(dataX, dataY)

MatTools.SaveXYData("WrightFisherWaitingTimeMoveU", dataX, dataY)

