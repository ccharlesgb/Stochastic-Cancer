# -*- coding: utf-8 -*-
"""
Created on Tue Feb 03 13:13:05 2015

@author: Connor
"""

import matplotlib.pyplot as plt
import Gillespie
import CarcinogenParam

myParam = CarcinogenParam.CarcinogenParam()
myGillespie = Gillespie.Gillespie()

myGillespie.params = myParam
myParam.Hook(myGillespie)

minr1 = 0.3
maxr1 = 3.0

DPC = 10
SDP = 500

dataX = []
dataY = []

for i in range(0, DPC):
    myParam.r1 = ((float(i) / (DPC - 1)) * (maxr1 - minr1)) + minr1
    dataX.append(myParam.r1)
    print(myParam.r1)
    fixCount = 0
    
    for i in range(0,SDP):
        myGillespie.Simulate()
        
        if myParam.n2 == myParam.N:
            fixCount += 1
    
    fixProb = float(fixCount) / SDP 
    dataY.append(fixProb)

plt.plot(dataX, dataY, linewidth=1.0, label="X2(t)")



