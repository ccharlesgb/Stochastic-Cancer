# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 13:52:16 2014

@author: Connor
"""

import SimTools
import matplotlib.pyplot as plt

import MatTools

mySim = SimTools.Gillespie(10)
mySim.in0 = 100

mySim.r0 = 1.0
mySim.r1 = 1.1
mySim.r2 = 0.9

mySim.u1 = 0.1
mySim.u2 = 0.1

mySim.timeLimit = 1000

mySim.populationHistory = 1

mySim.Simulate()

plt.plot(mySim.tHist, mySim.n0Hist)
plt.plot(mySim.tHist, mySim.n1Hist)
plt.plot(mySim.tHist, mySim.n2Hist)

plt.show()

MatTools.SaveRunHistory("BasicRunHist", mySim)

