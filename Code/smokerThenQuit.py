# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import SimTools
import time
import math

#Initialize the Gillespie simulator with 10 cells
mySim = SimTools.Gillespie(10)

mySim.in0 = 100
mySim.timeLimit = 100000.0
mySim.r1 = 1.0
mySim.r2 = 1.0

minSmokeTime = 0.0
maxSmokeTime = 1000.0

smokeTime = 0.0

smoke_u1 = 0.001
smoke_u2 = 0.1

quit_u1 = 0.0001
quit_u2 = 0.01

#Pre Simulate callback (called every frame before a timestep)
def increaseMutation(sim):
    global curPoint
    global dataPointCount
    global smokeTime

    smokeTime = (maxSmokeTime - minSmokeTime) * (float(curPoint) / (dataPointCount - 1)) + minSmokeTime
    
    if sim.curTime < smokeTime: #Still smoking
        sim.u1 = smoke_u1
        sim.u2 = smoke_u2
    else:
        sim.u1 = quit_u1
        sim.u2 = quit_u2

mySim.preSim = increaseMutation #IMPORTANT assign the callback (called in the class sim loop)
    
#Sweep the parameter r1 from 0.2 to 3.0 and run many simulations per data point
#Gets an idea on how likely cancer fixation is to occur for this parameter
simsPerDataPoint = 500
dataPointCount = 15

fixError = 1/math.sqrt(float(simsPerDataPoint))

#Initialize the array with default values
dataPointsX = []
dataPointsY = []
for i in range(0, dataPointCount):
    dataPointsX.insert(i,0.0)
    dataPointsY.insert(i,0.0)

maxSmokeTime = 2500

for curPoint in range(0, dataPointCount):
    startTime = time.clock() #Algorithm benchmarking
    
    fixationCounts = 0
    totFixTime = 0.0
    print("Current Data Point = {0}/{1} ({2}%)".format(curPoint + 1, dataPointCount, 100.0 * float(curPoint+1.0)/dataPointCount))
    #Perform many simulations to get an accurate probability of the fixation probability
    for sim in range(0, simsPerDataPoint):
        mySim.Simulate()
            
        totFixTime += mySim.curTime
        if mySim.Fixated() == 0:
            print("WARNING: DID NOT FIXATE")
            
    #Once the loop is done get the fraction of fixations for this r1
    dataPointsY[curPoint] = float(totFixTime) / float(simsPerDataPoint)
    if curPoint == 0:
        dataPointsX[curPoint] = smokeTime
        #maxSmokeTime = dataPointsY[0]
    else:
        dataPointsX[curPoint] = smokeTime

    print("Complete (Took {:.{s}f} seconds)".format(time.clock() - startTime, s=2))

#Dump data to console
for i in range(0, dataPointCount):
    print("radTime: {0} Fixation: {1}".format(dataPointsX[i],dataPointsY[i]))

#Create graph of data
plot_pulsed = plt.errorbar(dataPointsX, dataPointsY, yerr = fixError, label = "Smoking")
plt.xlabel("Smoke Time: ")
plt.ylabel("Type 2 Fixation Time")
plt.show()



        
