# -*- coding: utf-8 -*-
"""
Created on Thu Oct 23 11:41:06 2014

@author: Jonny
"""

import matplotlib.pyplot as plt
import TwoSpecies
import SimTools
import time
import math
import anfixtime

#Initialize the Gillespie simulator with 10 cells
N=50
mySim = TwoSpecies.Gillespie(N)
mySim.timeLimit = 1000
dataPointCount = N
mySim.j = int(N/2)
mySim.r0=1.0
mySim.r1=1.0


pulseOff = 0.01
pulseOn = 0.1
pulseWidth = 2.0
maxPulseWidth = 10.0

pulseWavelength = 10.0

#Pre Simulate callback (called every frame before a timestep)
def pulse_r1(sim):
    global curPoint
    global dataPointCount
    global pulseWidth

    pulseWidth = (float(curPoint) / (dataPointCount - 1)) * maxPulseWidth
    sim.r1 = SimTools.PulseWave(sim.curTime, pulseOn - pulseOff, pulseWidth, pulseWavelength) + pulseOff

#mySim.preSim = pulse_r1 #IMPORTANT assign the callback (called in the class sim loop)


def avgpulse_r1(sim):
    global curPoint
    global dataPointCount
    global pulseWidth
    pulseWidth = (float(curPoint) / (dataPointCount - 1)) * maxPulseWidth
    sim.r1 = (pulseOn - pulseOff) * (pulseWidth / pulseWavelength) + pulseOff
    
    
#Sweep the parameter r1 from 0.2 to 3.0 and run many simulations per data point
#Gets an idea on how likely cancer fixation is to occur for this parameter
simsPerDataPoint = 1000

#Initialize the array with default values
dataPointsX = []
dataPointsY = []
for i in range(0, dataPointCount):
    dataPointsX.insert(i,0.0)
    dataPointsY.insert(i,0.0)



for curPoint in range(0, dataPointCount):
    startTime = time.clock() #Algorithm benchmarking
    
    fixationTime = 0
    mySim.ij = int(float(curPoint)/(dataPointCount-1) * mySim.N)
    print("Current Data Point = {0}/{1} ({2}%)".format(curPoint + 1, dataPointCount, 100.0 * float(curPoint+1.0)/dataPointCount))
    #Perform many simulations to get an accurate probability of the fixation probability
    for sim in range(0, simsPerDataPoint):
        mySim.Simulate()
            
        if mySim.j == mySim.N or mySim.j==0: #The simulation ended with fixation
            fixationTime += mySim.curTime
    #Once the loop is done get the fraction of fixations for this r1
    dataPointsX[curPoint] = mySim.ij
    dataPointsY[curPoint] = float(fixationTime) / float(simsPerDataPoint)

    print("Complete (Took {:.{s}f} seconds)".format(time.clock() - startTime, s=2))

#Dump data to console
for i in range(0, dataPointCount):
    print("radTime: {0} Fixation: {1}".format(dataPointsX[i],dataPointsY[i]))

#Create graph of data
plot_pulsed = plt.errorbar(dataPointsX, dataPointsY, math.pow(simsPerDataPoint,0.5), label = "Pulsed")
plt.xlabel("r1 Pulse Time: ")
plt.ylabel("Fixation Time")
plt.show()

#Do the filename with all the parameters of the simulation   
filename = "TwoSpeciesPulsed_sim_N={0}_r0={1}_r1={2}_SPDP={3}".format(mySim.N, mySim.r0, mySim.r1, simsPerDataPoint)

#Save the data to a file
SimTools.SaveXYToFile(filename, dataPointsX, dataPointsY)

#Do the average of the pulse and see how it compares
#Initialize the array with default values
'''
#mySim.preSim = avgpulse_r1

dataPointsX = []
dataPointsY = []
for i in range(0, dataPointCount):
    dataPointsX.insert(i,0.0)
    dataPointsY.insert(i,0.0)

for curPoint in range(0, dataPointCount):
    startTime = time.clock() #Algorithm benchmarking
    
    fixationTime = 0
    
    print("Current Data Point = {0}/{1} ({2}%)".format(curPoint + 1, dataPointCount, 100.0 * float(curPoint+1.0)/dataPointCount))
    #Perform many simulations to get an accurate probability of the fixation probability
    for sim in range(0, simsPerDataPoint):
        mySim.Simulate()
            
        if mySim.Fixated(): #The simulation ended with fixation
            fixationTime += mySim.curTime
    #Once the loop is done get the fraction of fixations for this r1
    dataPointsX[curPoint] = pulseWidth
    dataPointsY[curPoint] = float(fixationTime) / float(simsPerDataPoint)

    print("Complete (Took {:.{s}f} seconds)".format(time.clock() - startTime, s=2))

#Dump data to console
for i in range(0, dataPointCount):
    print("radTime: {0} Fixation: {1}".format(dataPointsX[i],dataPointsY[i]))
'''
plot_averaged = plt.errorbar(dataPointsX, dataPointsY, math.pow(simsPerDataPoint,0.5), label = "Averaged")

plt.legend()

plt.show()

'''
#theoretical non-conditional fixation time (for averaged p)
dataPointsX = []
dataPointsY = []
for i in range(0,mySim.N):
    dataPointsX.append(i,i)
    dataPointsY.append(i,anfixtime.GetFixTimeJ(sim,i))

plt.plot(dataPointsX, dataPointsY, label = "Theoretical")'''