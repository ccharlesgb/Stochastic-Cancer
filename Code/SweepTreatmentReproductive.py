# -*- coding: utf-8 -*-
"""
Created on Thu Feb 05 13:58:42 2015

@author: Connor
"""

import GillespieTauLeap
import StemCellParam2
import matplotlib.pyplot as plt

myGill = GillespieTauLeap.Gillespie()

myParam = StemCellParam2.StemCellParam()

myHist = StemCellParam2.StemCellHist()

myGill.Hook(myParam) #Hook into sim parameters
myGill.SetHistory(myHist)

def EnableTreatment(time, params):
    if time > params.treatStart and time < params.treatStart + params.treatDur:
        params.rm = params.rm_ON
    else:
        params.rm = 0.0115
        
    if params.n[1] == 0:
        params.cureTime = time
        params.n[0] = 0
        params.n[1] = 0

myGill.preSim = EnableTreatment

myParam.n0[0] = 2e3
myParam.n0[1] = 5

myParam.rn = 0.005
myParam.rm = 0.0115

myParam.cn = 0.75e-3
myParam.cm = 0.38e-3
#myParam.cm = 0.75e-3

myParam.dn0 = 0.002
myParam.dm0 = 0.002

myParam.rn_ON = 0.002
myParam.rm_ON = 0.002

myParam.treatStart = 2000
myParam.treatDur = 1e100

myParam.cureTime = myGill.timeLimit

minMult = 0.01
maxMult = 0.3

dataX = []
dataY = []

pointCount = 8
simCount = 20

myGill.timeLimit = 1e10
myGill.tau = 2.0

for i in range(0, pointCount):
    cureTimeTot = 0.0
    print("Data Point: ", i)
    mult = float(maxMult - minMult) * float(i) / (pointCount - 1.0) + minMult
    print(mult)
    myParam.rm_ON = 0.0115 * mult
    for sim in range(0, simCount):
        myParam.cureTime = myGill.timeLimit  
        myHist.ClearFrames()
        myGill.Simulate()
        cureTimeTot += (myParam.cureTime - myParam.treatStart)
            
    dataX.append(mult)
    dataY.append(float(cureTimeTot) / simCount)

plt.plot(dataX, dataY, marker = 'o')
plt.yscale("log")
plt.xlabel("Reproductive multiplier")
plt.ylabel("Cure Time")