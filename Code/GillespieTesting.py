# -*- coding: utf-8 -*-
"""
Created on Tue Feb 03 13:13:05 2015

@author: Connor
"""

import matplotlib.pyplot as plt
import GillespieTauLeap
import CarcinogenParam

myParam = CarcinogenParam.CarcinogenParam()
myGillespie = GillespieTauLeap.Gillespie()

myGillespie.params = myParam
myParam.Hook(myGillespie)

myHist = CarcinogenParam.CarcinogenHist()

myGillespie.history = myHist

minr1 = 0.3
maxr1 = 3.0

DPC = 10
SDP = 500

dataX = []
dataY = []

myGillespie.timeLimit = 1000

myParam.addRate = 0.0

myParam.n0[0] = 1000
myParam.n0[1] = 0
myParam.n0[2] = 0

myGillespie.tau = 0.1

myParam.c0 = 0.003
myParam.c1 = 0.002
myParam.c2 = 0.001

myParam.u1 = 0.01
myParam.u2 = 0.01

myGillespie.Simulate()

'''
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
'''

plt.plot(myHist.tHist, myHist.n0Hist)
plt.plot(myHist.tHist, myHist.n1Hist)
plt.plot(myHist.tHist, myHist.n2Hist)
