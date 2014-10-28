# -*- coding: utf-8 -*-
"""
Created on Thu Oct 23 16:46:30 2014

@author: Jonny
"""

#theoretical non-conditional fixation time (for averaged p)
import matplotlib.pyplot as plt
import TwoSpecies
import anfixtime

#Initialize the Gillespie simulator with 10 cells
N=10
mySim = TwoSpecies.Gillespie(N)
mySim.timeLimit = 1000000
dataPointCount = 25
mySim.j = int(N/2)
mySim.r0=1.1
mySim.r1=1.5

dataPointsX = []
dataPointsY = []
for i in range(1,mySim.N):
    dataPointsX.append(i)
    dataPointsY.append(anfixtime.GetFixTimeJ(mySim,i))

plt.plot(dataPointsX, dataPointsY, label = "Theoretical")