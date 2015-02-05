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
#set some paramters
myParam.in0 = 2e3
myParam.im0 = 1

myParam.rn = 0.005
myParam.rm = 0.0115

myParam.cn = 0.75e-3
#myParam.cm = 0.38e-3
myParam.cn = 0.75e-3

myParam.dn0 = 0.002
myParam.dm0 = 0.002

myParam.dm0_ON = 0.002
myParam.dm0_ON = 0.002

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