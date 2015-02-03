# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 16:10:45 2014
@author: Jonny
"""

import systematic_transitions
import SimTools
import numpy as np
import math
import matplotlib.pyplot as plt
import MatTools
import DetermenisticModel

#y=np.zeros((10,10))
#print(y[0][0])
#initialise simulation

#general simulation parameters
mySim = SimTools.Gillespie(10)

mySim.in0 = 10
mySim.timeLimit = 10000.0
mySim.u1 = 0.1
mySim.u2 = 0.1
mySim.r0 = 1.0
mySim.r1 = 1.0
mySim.r2 = 1.0

mySim.ResetSim()

deterSim = DetermenisticModel.Determenistic()

dataPointCount = 5
spd = 100

min_r1 = 0.5
max_r1 = 1.5

min_u1 = 0.01
max_u1 = 0.1

min_N = 10
max_N = 50

dataX = []
dataFix = []
dataTheory = []

mySim.ResetSim()

for curPoint in range(0,dataPointCount): #simulate some shit
    #mySim.r1 = float(max_r1 - min_r1) * float(curPoint) / (dataPointCount - 1.0) + min_r1
    #mySim.u1 = float(max_u1 - min_u1) * float(curPoint) / (dataPointCount - 1.0) + min_u1
    mySim.in0 = int(max_N - float(max_N - min_N) * float(curPoint) / (dataPointCount - 1.0))
    dataX.append(mySim.in0)
    print(mySim.in0)
    print("Data Point {0}/{1} - {2}%".format(curPoint + 1, dataPointCount, 100*(curPoint + 1)/dataPointCount))
    
    fixTimeTerm = 0.0
    for i in range(0,spd):
        mySim.Simulate()
        if mySim.Fixated():
            fixTimeTerm += mySim.curTime
        else:
            fixTimeTerm += mySim.timeLimit
            print("DID NOT FIXATE")
    dataFix.append(fixTimeTerm/spd)
    
    #Do theory
    mySim.ResetSim()    
    deterSim.Reset()
    deterSim.Integrate(mySim)
    
    result = deterSim.curT
    dataTheory.append(result)




#MatTools.SaveXYData("Jnumerical_fixation_time_varying_r1_02_to_30", r1, fixTime_num, xLabel = "r1", yLabel = "Fixation Time", sim = mySim)
#MatTools.SaveXYData("Jsimulation_fixation_time_varying_r1_02_30", r1, fixTime_num, yError = errors ,xLabel = "r1", yLabel = "Fixation Time", sim = mySim)

plt.plot(dataX,dataFix, '^')
plt.plot(dataX,dataTheory)
plt.xlabel("N")
plt.ylabel("Fixation Time")
plt.show()

file_name = "ThreeSpeciesDeterm_sweepN_N_{0}_SDP_{1}".format(mySim.N, spd)

data2 = dict()
data2["num_fixtime_N"] = dataTheory

MatTools.SaveXYData(file_name, dataX, dataFix, xLabel = "N", yLabel = "fixtime_N", sim = mySim, otherDict = data2)