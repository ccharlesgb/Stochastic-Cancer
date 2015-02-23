# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 18:53:15 2015

@author: Jonny
"""

import wright_fisher
import matplotlib.pyplot as plt

cellTypes = 60
population = 1000000

myWF = wright_fisher.wright_fisher()
myHist = wright_fisher.wf_hist(cellTypes)

myWF.s = 0.1

myWF.cellTypes = cellTypes
myWF.popSize = population

myWF.history = myHist

myWF.popSize = population

myWF.stepLimit = 3000

myWF.Simulate()

for i in range(0, cellTypes):
   plt.plot(myHist.stepHist, myHist.histArray[i])
