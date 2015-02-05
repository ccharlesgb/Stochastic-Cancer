# -*- coding: utf-8 -*-
"""
Created on Thu Feb  5 16:09:36 2015

@author: Jonny
"""


import Gillespie
import StemCellParam2
import matplotlib.pyplot as plt

myGill = Gillespie.Gillespie()

myParam = StemCellParam2.StemCellParam()

myHist = StemCellParam2.StemCellHist()

myGill.Hook(myParam) #Hook into sim parameters
myGill.SetHistory(myHist)

DPC = 1
print("GET READY BEFORE")
myGill.Simulate()

plt.subplot(2, 1, 1)
plt.plot(myHist.tHist, myHist.n0Hist)
plt.plot(myHist.tHist, myHist.m0Hist)
plt.ylabel('Number of Cells - stem')

plt.subplot(2, 1, 2)
plt.plot(myHist.tHist, myHist.n1Hist)
plt.plot(myHist.tHist, myHist.m1Hist)
plt.xlabel('time (Days)')
plt.ylabel('Number of Cells - differentiated')

plt.show()