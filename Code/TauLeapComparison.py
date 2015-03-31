# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 16:23:16 2015

@author: Connor
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Mar 09 15:09:08 2015

@author: Connor
"""

import matplotlib.pyplot as plt
import MatTools
import numpy as np
import math
import time
import DeterministicTauLeap

import TauLeap
import TauLeapParam
import wright_fisher

cellTypes = 3
timeLimit = 1e4

myParam = TauLeapParam.Params(cellTypes)

#Create the simulator
myGillespie = TauLeap.Sim(cellTypes)
myGillespie.params = myParam
myParam.Hook(myGillespie)
#Set History
myHist = TauLeapParam.Hist(cellTypes)
myGillespie.history = myHist

myGillespie.timeLimit = timeLimit
myGillespie.RECORD_TAU_INFO = 0
myGillespie.n_c = 10
myGillespie.stopAtAppear = 0
myGillespie.epsilon = 0.2

#Param Stuff
myParam.n0[0] = 1e4
s = 0.01
for i in range(0,cellTypes):
    myParam.r[i] = math.pow(1.0 + s, i)
    myParam.u[i] = 1e-7
myParam.USE_D = True
myParam.d = 100

#Deter Stuff
deterSim = DeterministicTauLeap.Sim(cellTypes)
deterSim.params = myParam
myHist2 = TauLeapParam.Hist(cellTypes)
deterSim.history = myHist2

deterSim.timeStep = 0.01
deterSim.timeLimit = timeLimit

myWF = wright_fisher.wright_fisher()

myGillespie.Simulate()  
deterSim.Integrate()

plt.figure()
for i in range(0, cellTypes):
   plt.plot(myHist.tHist, myHist.histArray[i])
   plt.plot(myHist2.tHist, myHist2.histArray[i], ':')
   plt.yscale("log")
