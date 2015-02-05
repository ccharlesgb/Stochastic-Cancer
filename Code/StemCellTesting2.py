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
        params.dm0 = 0.0125
        params.dn0 = 0.003
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

simCount = 5
for i in range(0, simCount):
    myHist.ClearFrames()
    myGill.Simulate()
    
    print ("DONE")
    
    if myParam.m0 == 0:    
    
    plt.title("Stem Cells")
    plt.plot(myHist.tHist, myHist.n0Hist)
    plt.plot(myHist.tHist, myHist.m0Hist, '--')
#plt.yscale("log")