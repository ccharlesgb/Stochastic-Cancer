# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 11:26:29 2015

@author: Jonny
"""

import matplotlib.pyplot as plt
import GillespieTauLeap
import CarcinogenParam
import gillespie

myParam = CarcinogenParam.CarcinogenParam()
tauGillespie = GillespieTauLeap.Gillespie()
exactGillespie = gillespie.Gillespie()

tauGillespie.params = myParam
exactGillespie.params = myParam
tauGillespie.Hook(myParam)
exactGillespie.Hook(myParam)

tauHist = CarcinogenParam.CarcinogenHist()
exactHist = CarcinogenParam.CarcinogenHist()

tauGillespie.history = tauHist
exactGillespie.history = exactHist

minr1 = 0.3
maxr1 = 3.0

DPC = 10
SDP = 500

dataX_tau = []
dataY_tau = []
dataX_exact = []
dataY_exact = []

tauGillespie.timeLimit = 100

myParam.addRate = 0.0

myParam.n0[0] = 20
myParam.n0[1] = 0
myParam.n0[2] = 0

tauGillespie.tau = 0.1

myParam.c0 = 0.003
myParam.c1 = 0.002
myParam.c2 = 0.001

myParam.u1 = 0.1
myParam.u2 = 0.1

#simulate using tau leaping
for i in range(0, DPC):
    myParam.r1 = ((float(i) / (DPC - 1)) * (maxr1 - minr1)) + minr1
    dataX_tau.append(myParam.r1)
    print(myParam.r1)
    fixCount = 0
    
    for i in range(0,SDP):
        tauGillespie.Simulate()
        
        if myParam.n[2] == myParam.N:
            fixCount += 1
    
    fixProb = float(fixCount) / SDP 
    dataY_tau.append(fixProb)

#simulate using SSA
for i in range(0, DPC):
    myParam.r1 = ((float(i) / (DPC - 1)) * (maxr1 - minr1)) + minr1
    dataX_exact.append(myParam.r1)
    print(myParam.r1)
    fixCount = 0
    
    for i in range(0,SDP):
        tauGillespie.Simulate()
        
        if myParam.n[2] == myParam.N:
            fixCount += 1
    
    fixProb = float(fixCount) / SDP 
    dataY_exact.append(fixProb)

plt.plot(dataX_tau, dataY_tau, marker = 'o', label="tau leaping")
plt.plot(dataX_exact, dataY_exact, linewidth=1.0, label="tau leaping")
plt.legend()
plt.show()