# -*- coding: utf-8 -*-
"""
Created on Tue Feb 10 13:53:04 2015

@author: Connor
"""

import NGillespie
import NTypeParam
import matplotlib.pyplot as plt

tCount = 10

myGill = NGillespie.Gillespie()
myParam = NTypeParam.NParam(tCount)

myHist = NTypeParam.NHist(tCount)
myGill.history = myHist

myParam.ini[0] = 100

myGill.timeLimit = 10000

myGill.params = myParam

myGill.Simulate()

for i in range(0, tCount):
   plt.plot(myHist.tHist, myHist.histArray[i])
