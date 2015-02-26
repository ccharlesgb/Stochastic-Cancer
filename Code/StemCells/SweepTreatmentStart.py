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
    if time > myParam.treatStart and time < myParam.treatStart + myParam.treatDur:
        params.dm0 = params.dm0_ON
        #params.dn0 = params.dn0_ON
    else:
        params.dm0 = 0.002
        #params.dn0 = 0.002
        
    if params.n[1] == 0:
        params.cureTime = time
        params.n[0] = 0
        params.n[1] = 0

myGill.preSim = EnableTreatment

myGill.timeLimit = 1e10

myParam.n0[0] = 2e3
myParam.n0[1] = 5

myParam.rn = 0.005
myParam.rm = 0.0115

myParam.cn = 0.75e-3
myParam.cm = 0.38e-3
#myParam.cm = 0.75e-3

myParam.dn0 = 0.002
myParam.dm0 = 0.002

myParam.treatStart = 2000
myParam.treatDur = 1e10

myParam.dn0_ON = 0.002 * 7.5
myParam.dm0_ON = 0.002 * 7.5

myParam.cureTime = myGill.timeLimit

minTime = 1000
maxTime = 4000

dataX = []
dataY = []

pointCount = 8
simCount = 200

myGill.timeLimit = 1e5
myGill.tau = 3.0

for i in range(0, pointCount):
    cureTimeTot = 0.0
    print("Data Point: ", i)
    startTime = float(maxTime - minTime) * float(i) / (pointCount - 1.0) + minTime
    print(startTime)  
    
    for sim in range(0, simCount):
        myParam.treatStart = startTime 
        myParam.cureTime = myGill.timeLimit  
        myHist.ClearFrames()
        myGill.Simulate()
        cureTimeTot += (myParam.cureTime - myParam.treatStart)
            
    dataX.append(myParam.treatStart)
    dataY.append(float(cureTimeTot) / simCount)

plt.plot(dataX, dataY, marker = 'o')
plt.xlabel("Treatment Start Time")
plt.ylabel("Cure Time")