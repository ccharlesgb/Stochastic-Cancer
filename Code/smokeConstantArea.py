# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import SimTools
import time
import math

#Initialize the Gillespie simulator with 10 cells
mySim = SimTools.Gillespie(10)

mySim.in0 = 10
mySim.timeLimit = 10000.0
mySim.r1 = 1.0
mySim.r2 = 1.0

minSmokeTime = 0.0
maxSmokeTime = 250

smokeArea_u2 = 2.0
smokeArea_u1 = 2.0

quit_u1 = 0.002
quit_u2 = 0.002

smoke_u1 = quit_u1 * 10.0
smoke_u2 = quit_u2 * 10.0

min_mut_factor = 2.0
max_mut_factor = 20.0

#Pre Simulate callback (called every frame before a timestep)
def increaseMutation(sim):
    global curPoint
    global dataPointCount
    global smokeTime
  
    if sim.curTime < smokeTime: #Still smoking
        sim.u1 = smoke_u1
        sim.u2 = smoke_u2
    else:
        sim.u1 = quit_u1
        sim.u2 = quit_u2

mySim.preSim = increaseMutation #IMPORTANT assign the callback (called in the class sim loop)
    
#Sweep the parameter r1 from 0.2 to 3.0 and run many simulations per data point
#Gets an idea on how likely cancer fixation is to occur for this parameter
simsPerDataPoint = 5000
dataPointCount = 5

fixError = []

#Initialize the array with default values
dataPointsX = []
dataPointsY = []

dataST = []
dataMU = []

MUS = []
for i in range(0, dataPointCount):
    val = (float(i) / (dataPointCount - 1))

for curPoint in range(0, dataPointCount):
    startTime = time.clock() #Algorithm benchmarking
    
    fixationCounts = 0
    totFixTime = 0.0
    print("Current Data Point = {0}/{1} ({2}%)".format(curPoint + 1, dataPointCount, 100.0 * float(curPoint+1.0)/dataPointCount))
    
    #mut_factor = (max_mut_factor - min_mut_factor) * ((dataPointCount) / float(curPoint + 1.0)) + min_mut_factor
    mut_factor = (max_mut_factor - min_mut_factor) * ((1.0) / float(curPoint + 1.0)) + min_mut_factor
    smoke_u1 = quit_u1 * mut_factor
    smokeTime = smokeArea_u1 / (smoke_u1 - quit_u1)

    smoke_u2 = smoke_u1
    
    dataST.append(smokeTime)
    dataMU.append(smoke_u1)
    
    print("MUT FACTOR: {0}".format(smoke_u1))
    print("T SMOKE: {0}".format(smokeTime))
    
    #Perform many simulations to get an accurate probability of the fixation probability
    for sim in range(0, simsPerDataPoint):
        mySim.Simulate()
        
        totFixTime += mySim.curTime
        if mySim.Fixated() == 0:
            print("WARNING NO FIX")
        fixationCounts += mySim.Fixated()
    
    #Once the loop is done get the fraction of fixations for this r1
    #dataPointsY.append(float(fixationCounts) / float(simsPerDataPoint))

    fixTimeAvg = float(totFixTime) / float(simsPerDataPoint)
    
    fixError.append(fixTimeAvg/math.sqrt(float(simsPerDataPoint)))    
    
    dataPointsY.append(fixTimeAvg)
    dataPointsX.append(smokeTime)
    #maxSmokeTime = dataPointsY[0]

    print("Complete (Took {:.{s}f} seconds)".format(time.clock() - startTime, s=2))

#Dump data to console
for i in range(0, len(dataPointsX)):
    print("radTime: {0} Fixation: {1}".format(dataPointsX[i],dataPointsY[i]))

plt.figure()
plt.subplot(2,1,1)

plt.plot(dataST, dataMU, 'o-')
plt.ylabel("Mutation Amount")
plt.xlabel("Smoke Time")

#Create graph of data+
plt.subplot(2,1,2)
plot_pulsed = plt.errorbar(dataPointsX, dataPointsY, yerr = fixError, label = "Smoking")
plt.xlabel("Smoke Time: ")
plt.ylabel("Type 2 Fixation Time")
plt.show()



        
