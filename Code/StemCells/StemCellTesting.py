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

myGill.timeLimit = 50

myParam.in0 = 10
myParam.in1 = 100
myParam.im0 = 5
myParam.im1 = 0

myParam.rn = 0.05e2
myParam.rm = 0.115e2

myParam.cn = 0.75e1
#myParam.cm = 0.38e-4
myParam.cm = 0.75e1

myParam.dn0 = 0.002
myParam.dn1 = 0.213
myParam.dm0 = 0.002
myParam.dm1 = 0.213

myParam.an = 1.065e3
myParam.am = 1.065e3

myGill.Simulate()

print ("DONE")

plt.subplot(211)
plt.title("Stem Cells")
plt.plot(myHist.tHist, myHist.n0Hist)
plt.plot(myHist.tHist, myHist.m0Hist, '--')
#plt.yscale("log")

plt.subplot(212)
plt.title("Diff Cells")
plt.plot(myHist.tHist, myHist.n1Hist)
plt.plot(myHist.tHist, myHist.m1Hist, '--')
#plt.yscale("log")
