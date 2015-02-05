# -*- coding: utf-8 -*-
"""
Created on Thu Feb 05 13:58:42 2015

@author: Connor
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Feb 03 14:15:29 2015

@author: Connor
"""

import Gillespie
import StemCellParam2
import matplotlib.pyplot as plt

myGill = Gillespie.Gillespie()

myParam = StemCellParam2.StemCellParam()

myHist = StemCellParam2.StemCellHist()

myGill.Hook(myParam) #Hook into sim parameters
myGill.SetHistory(myHist)

def EnableTreatment(time, params):
    treatStart = 2000.0
    treatDur = 2000
    if time > treatStart and time < treatStart + treatDur:
        params.dm0 = params.dm0_ON
        params.dn0 = params.dm0_ON
    else:
        params.dm0 = 0.002
        params.dn0 = 0.002

myGill.preSim = EnableTreatment

myGill.timeLimit = 5000

myParam.in0 = 2e3
myParam.im0 = 5

myParam.rn = 0.005
myParam.rm = 0.0115

myParam.cn = 0.75e-3
#myParam.cm = 0.38e-3
myParam.cn = 0.75e-3

myParam.dn0 = 0.002
myParam.dm0 = 0.002

myParam.dn0_ON = 0.002
myParam.dm0_ON = 0.002

minMult = 5
maxMult = 10

dataX = []
dataY = []

pointCount = 30
simCount = 10

for i in range(0, pointCount):
    cureCount = 0
    print("Data Point: ", i)
    mult = float(maxMult - minMult) * float(i) / (pointCount - 1.0) + minMult
    print(mult)
    for sim in range(0, simCount):
        myParam.dm0_ON = mult * 0.002   
        
        myHist.ClearFrames()
        myGill.Simulate()
        
        if myParam.m0 == 0:
            cureCount += 1
            
    dataX.append(myParam.dm0_ON)
    dataY.append(float(cureCount) / simCount)

plt.plot(dataX, dataY, marker = 'o')
plt.xlabel("Treatment death rate")
plt.ylabel("Cure Probability")