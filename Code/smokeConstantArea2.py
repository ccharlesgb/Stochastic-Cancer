# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import anfixtime
import TwoSpecies
import time
import math

import MatTools

#Initialize the Gillespie simulator with 10 cells
mySim = TwoSpecies.Gillespie(10)

mySim.ij = 0
mySim.N = 10
mySim.timeLimit = 10000.0
mySim.r0 = 1.0
mySim.r1 = 1.0

minSmokeTime = 0.0
maxSmokeTime = 250

smokeArea_u1 = 1.0

quit_u1 = 0.01

smoke_u1 = quit_u1 * 10.0

min_mut_factor = 2.0
max_mut_factor = 20.0

#Pre Simulate callback (called every frame before a timestep)
def increaseMutation(sim):
    global curPoint
    global dataPointCount
    global smokeTime
  
    if sim.curTime < smokeTime: #Still smoking
        sim.u1 = smoke_u1
    else:
        sim.u1 = quit_u1

mySim.preSim = increaseMutation #IMPORTANT assign the callback (called in the class sim loop)
    
#Sweep the parameter r1 from 0.2 to 3.0 and run many simulations per data point
#Gets an idea on how likely cancer fixation is to occur for this parameter
simsPerDataPoint = 100000
dataPointCount = 10

fixError = []

#Initialize the array with default values
dataPointsX = []
dataPointsY = []

dataST = []
dataMU = []


for curPoint in range(0, dataPointCount):
    startTime = time.clock() #Algorithm benchmarking
    
    fixationCounts = 0
    totFixTime = 0.0
    print("Current Data Point = {0}/{1} ({2}%)".format(curPoint + 1, dataPointCount, 100.0 * float(curPoint+1.0)/dataPointCount))
    
    #mut_factor = (max_mut_factor - min_mut_factor) * ((dataPointCount) / float(curPoint + 1.0)) + min_mut_factor
    mut_factor = (max_mut_factor - min_mut_factor) * ((1.0) / float(curPoint + 1.0)) + min_mut_factor
    smoke_u1 = quit_u1 * mut_factor
    smokeTime = smokeArea_u1 / (smoke_u1 - quit_u1)
    
    dataST.append(smokeTime)
    dataMU.append(smoke_u1)
    
    print("MUT FACTOR: {0}".format(smoke_u1))
    print("T SMOKE: {0}".format(smokeTime))
    
    print("CHECK DOSE: {0}".format((smoke_u1 - quit_u1) * smokeTime))
    
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
plt.ylabel("u1 / u2")
plt.xlabel("t_mut")

MatTools.SaveXYData("MutIncrease_Params", dataST, dataMU, xLabel = "SmokeTime", yLabel = "u1 & u2", sim = mySim)

#Create graph of data+
plt.subplot(2, 1, 2)
plot_pulsed = plt.errorbar(dataPointsX, dataPointsY, yerr = fixError)
plt.xlabel("t_mut")
plt.ylabel("Type-1 Fixation Time")
plt.show()

MatTools.SaveXYData("MutIncrease_FixTimes", dataPointsX, dataPointsY, xLabel = "SmokeTime", yLabel = "Fixation Time", sim = mySim)
