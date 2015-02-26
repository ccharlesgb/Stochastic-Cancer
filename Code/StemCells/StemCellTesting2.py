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

import GillespieTauLeap
import StemCellParam2
import matplotlib.pyplot as plt

myGill = GillespieTauLeap.Gillespie()

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

myParam.n0[0] = 2e3
myParam.n0[1] = 2

myParam.rn = 0.005
myParam.rm = 0.0115

myParam.cn = 0.75e-3
#myParam.cm = 0.38e-3
myParam.cn = 0.75e-3

myParam.dn0 = 0.002
myParam.dm0 = 0.002

myParam.dm0_ON = 0.002
myParam.dm0_ON = 0.002

minMult = 1.1
maxMult = 1.5

dataX = []
dataY = []


simCount = 10

for sim in range(0, simCount):
    mult = 7
    print("DONE")
    myParam.dm0_ON = mult * 0.002   
    
    myHist.ClearFrames()
    myGill.Simulate()
    
    print(myParam.n[1])

    plt.plot(myHist.tHist, myHist.n0Hist)
    plt.plot(myHist.tHist, myHist.m0Hist, '--')

