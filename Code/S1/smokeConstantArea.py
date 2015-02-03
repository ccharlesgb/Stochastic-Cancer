# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import SimTools
import time
import math
import MatTools

#Initialize the Gillespie simulator with 10 cells
mySim = SimTools.Gillespie(10)

mySim.in0 = 10
mySim.timeLimit = 200000.0
mySim.r1 = 1.0
mySim.r2 = 1.0

smokeArea_u = 5.0

quit_u = 0.001

smoke_u = quit_u * 10.0

min_mut_factor = 5.0
max_mut_factor = 500.0

#Pre Simulate callback (called every frame before a timestep)
def increaseMutation(sim):
    global curPoint
    global dataPointCount
    global smokeTime
  
    if sim.curTime < smokeTime: #Still smoking
        sim.u1 = smoke_u
        sim.u2 = smoke_u
    else:
        sim.u1 = quit_u
        sim.u2 = quit_u

mySim.preSim = increaseMutation #IMPORTANT assign the callback (called in the class sim loop)
    
#Sweep the parameter r1 from 0.2 to 3.0 and run many simulations per data point
#Gets an idea on how likely cancer fixation is to occur for this parameter
simsPerDataPoint = int(2e6)
dataPointCount = 6
fixError = []

#Initialize the array with default values
dataPointsX = []
dataPointsY = []

dataST = []
dataMU = []

MUS = []
for i in range(0, dataPointCount):
    val = (float(i) / (dataPointCount - 1))

minSmokeT = 10
maxSmokeT = 80

for curPoint in range(0, dataPointCount):
    startTime = time.clock() #Algorithm benchmarking
    
    fixationCounts = 0
    totFixTime = 0.0
    print("Current Data Point = {0}/{1} ({2}%)".format(curPoint + 1, dataPointCount, 100.0 * float(curPoint+1.0)/dataPointCount))
    
    smokeTime = (maxSmokeT - minSmokeT) * (curPoint / (dataPointCount - 1.0)) + minSmokeT
    smoke_u = smokeArea_u / smokeTime + quit_u
    
    '''
    #mut_factor = (max_mut_factor - min_mut_factor) * ((dataPointCount) / float(curPoint + 1.0)) + min_mut_factor
    mut_factor = (max_mut_factor - min_mut_factor) * ((1.0) / float(curPoint + 1.0)) + min_mut_factor
    smoke_u = quit_u * mut_factor
    smokeTime = smokeArea_u / (quit_u * (mut_factor - 1.0))
    '''
    dataST.append(smokeTime)
    dataMU.append(smoke_u)
    
    print("MUT FACTOR: {0}".format(smoke_u))
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
    
    if fixTimeAvg < smokeTime:
        print("WARNING FIX TIME < SMOKE TIME")
    
    dataPointsY.append(fixTimeAvg)
    dataPointsX.append(smokeTime)
    #maxSmokeTime = dataPointsY[0]

    print("Complete (Took {:.{s}f} seconds)".format(time.clock() - startTime, s=2))

#Dump data to console
for i in range(0, len(dataPointsX)):
    print("radTime: {0} Fixation: {1}".format(dataPointsX[i],dataPointsY[i]))


def PlotPulse(area, mut_factor, offset):
    dataX = []
    dataY = []
    
    mutTime = area / (offset * (mut_factor - 1.0))
    dataX.append(0.0)
    dataX.append(mutTime)
    dataX.append(mutTime)
    dataX.append(100.0)
    dataY.append(offset * mut_factor)
    dataY.append(offset * mut_factor)
    dataY.append(offset)
    dataY.append(offset)
    
    plt.plot(dataX, dataY)
    

plt.figure()
plt.subplot(2,1,1)

for i in range(len(dataST)):
    PlotPulse(smokeArea_u, dataMU[i] / quit_u, quit_u)

#plt.plot(dataST, dataMU, 'o-')
#plt.ylabel("u1 = u2 = u")
#plt.xlabel("Pulse Time")

#plt.plot()

#Create graph of data+
plt.subplot(2,1,2)
plot_pulsed = plt.errorbar(dataPointsX, dataPointsY, yerr = fixError, label = "Smoking")
plt.xlabel("Pulse Time")
plt.ylabel("Type 2 Fixation Time")
plt.show()

file_name = "PulseConstantArea_r1_{0}_r2_{1}_uquit_{2}_N_{3}_Dose_{4}_SDP_{5}".format(mySim.r1, mySim.r2, quit_u, mySim.N, smokeArea_u, simsPerDataPoint)
MatTools.SaveXYData("PARAM_" + file_name, dataST, dataMU, xLabel = "PulseTime", yLabel = "u_pulse")
MatTools.SaveXYData("DAT_" + file_name, dataPointsX, dataPointsY, xLabel = "PulseTime", yLabel = "FixTime")



        
