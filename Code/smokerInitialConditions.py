# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import SimTools
import time
import math

#Initialize the Gillespie simulator with 10 cells
mySim = SimTools.Gillespie(10)

mySim.in0 = 100
mySim.N = 100
mySim.timeLimit = 100.0
mySim.r1 = 1.0
mySim.r2 = 1.0
mySim.u1 = 0.01
mySim.u2 = 0.01
    
#Sweep the parameter r1 from 0.2 to 3.0 and run many simulations per data point
#Gets an idea on how likely cancer fixation is to occur for this parameter
simsPerDataPoint = 5000
dataPointCount = 10

fixError = 1/math.sqrt(float(simsPerDataPoint))

#Initialize the array with default values
dataPointsX = []
dataPointsY = []

minn1 = 0.0
maxn1 = float(mySim.N) * 1.0

for curPoint in range(0, dataPointCount):
    startTime = time.clock() #Algorithm benchmarking

    IN1 = (maxn1 - minn1) * float(curPoint) / (dataPointCount - 1) + minn1
    IN1 = int(IN1)
    print("IN1 = {0}".format(IN1))
    mySim.in1 = IN1
    mySim.in0 = mySim.N - IN1
    
    fixationCounts = 0
    totFixTime = 0.0
    print("Current Data Point = {0}/{1} ({2}%)".format(curPoint + 1, dataPointCount, 100.0 * float(curPoint+1.0)/dataPointCount))
    #Perform many simulations to get an accurate probability of the fixation probability
    for sim in range(0, simsPerDataPoint):
        mySim.Simulate()
        
        fixationCounts += mySim.Fixated()            
            
        totFixTime += mySim.curTime
        #if mySim.Fixated() == 0:
            #print("WARNING: DID NOT FIXATE")
    
    #Once the loop is done get the fraction of fixations for this r1
    #dataPointsY.append(float(totFixTime) / float(simsPerDataPoint))
    dataPointsY.append(float(fixationCounts) / simsPerDataPoint)
    dataPointsX.append(IN1)
    #maxSmokeTime = dataPointsY[0]
    
    print("Complete (Took {:.{s}f} seconds)".format(time.clock() - startTime, s=2))

#Dump data to console
for i in range(0, len(dataPointsX)):
    print("radTime: {0} Fixation: {1}".format(dataPointsX[i],dataPointsY[i]))

#Create graph of data
plt.figure()
plot_pulsed = plt.errorbar(dataPointsX, dataPointsY, yerr = fixError, label = "Smoking")
plt.xlabel("Smoke Time: ")
plt.ylabel("Type 2 Fixation Time")
plt.show()


        
