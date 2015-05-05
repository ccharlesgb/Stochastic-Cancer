# -*- coding: utf-8 -*-
"""
Created on Tue May 05 10:40:34 2015

@author: Connor
"""

"""
A simple example of an animated plot
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math

import wright_fisher
import TauLeapParam


def GetDistribution(t):
    dataY = []
    for i in range(0, cellTypes): 
        dataY.append(myHist.histArray[i][t])
    return dataY

def GetAppearanceTimes():
    for i in range(0, cellTypes):
        found = False     
        for t in range(len(myHist.histArray[i])):
            if found:
                continue
            count = myHist.histArray[i][t]
            if count > 0:
                avgAppear[i] += t
                found = True

cellTypes = 6
population = 1e9

myWF = wright_fisher.wright_fisher(cellTypes)
myHist = TauLeapParam.Hist(cellTypes)
myParam = TauLeapParam.Params(cellTypes)

myWF.stopAtAppear = 1
myWF.history = myHist
myWF.stepLimit = 20000
myWF.useApproxTheta = 0
myWF.params = myParam

myParam.n0[0] = population
myParam.d = 100
myParam.uNotConst = 0

s = 1e-2

myParam.SetCompoundFitness(1e-2)
myParam.SetUAll(1e-7)

myWF.Simulate()
    
avgAppear = [0.0] * cellTypes

vel_hist = []
avg_hist = []

#Plotting    
LINE_HAND = []

fig, ax = plt.subplots(figsize = [16,9],facecolor='white')
for i in range(0,cellTypes):
    line, = ax.semilogy(myHist.tHist[0:10], myHist.histArray[i][0:10], '-')
    LINE_HAND.append(line)

ax.set_xlim(0,myWF.curTime)
ax.set_ylim(1, 1e4)
ax.set_xlabel("t")
ax.set_ylabel("N_j")

#plt.legend(frameon=False, loc = 2)

def animate(t):
    for lnum,line in enumerate(LINE_HAND):
        line.set_xdata(myHist.tHist[0:t])
        line.set_ydata(myHist.histArray[lnum][0:t])
    return tuple(LINE_HAND)

#Init only required for blitting to give a clean slate.
def init():
    for lnum,line in enumerate(LINE_HAND):
        line.set_xdata(np.ma.array(myHist.tHist[0:1], mask=True))
        line.set_ydata(np.ma.array(myHist.histArray[i][0:1], mask=True))
    return tuple(LINE_HAND)

ani = animation.FuncAnimation(fig, animate, np.append([0]*10,np.append(np.arange(0, myWF.curTime,3),[myWF.curTime]*10)), init_func=init,
    interval=20, blit=True)
plt.show()