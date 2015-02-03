# -*- coding: utf-8 -*-
"""
Created on Tue Feb 03 14:15:29 2015

@author: Connor
"""

import Gillespie
import StemCellParam
import matplotlib.pyplot as plt

myGill = Gillespie.Gillespie()

myParam = StemCellParam.StemCellParam()

myHist = StemCellParam.StemCellHist()

myGill.Hook(myParam) #Hook into sim parameters
myGill.SetHistory(myHist)

myGill.Simulate()

print ("DONE")

plt.plot(myHist.tHist, myHist.n0Hist)
plt.plot(myHist.tHist, myHist.n1Hist, ':')
plt.plot(myHist.tHist, myHist.m0Hist, '--')
plt.plot(myHist.tHist, myHist.m1Hist)
