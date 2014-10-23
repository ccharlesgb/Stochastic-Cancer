# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import SimTools
import time
import math

#Initialize the Gillespie simulator with 10 cells
mySim = SimTools.Gillespie(100)
mySim.timeLimit = 1000.0
mySim.u1 = 0.1
mySim.u2 = 0.1
mySim.r1 = 0.75
mySim.populationHistory = 1

dataPointCount = 10

pulseAmp = mySim.r1
pulseWidth = 5.0
maxPulseWidth = 10.0

pulseWavelength = 10.0

dataPulseTime = []
dataPulsePulse = []

#Pre Simulate callback (called every frame before a timestep)
def pulse_r1(sim):
    global curPoint
    global dataPointCount
    global pulseWidth

    sim.r1 = SimTools.PulseWave(sim.curTime, pulseAmp, pulseWidth, pulseWavelength) + pulseAmp
    
    dataPulseTime.append(sim.curTime)
    dataPulsePulse.append(sim.r1)  
    
mySim.preSim = pulse_r1 #IMPORTANT assign the callback (called in the class sim loop)

#Run the simulation
mySim.Simulate()

plt.subplot(2,1,1)
plt.plot(dataPulseTime, dataPulsePulse)
plt.xlabel("time")
plt.ylabel("r1")

#Create graph of data
plt.subplot(2, 1, 2)
plot_pulsed = plt.plot(mySim.tHist, mySim.n0Hist, label = "n0", marker = 'o')
plot_pulsed = plt.plot(mySim.tHist, mySim.n1Hist, label = "n1", marker = '^')
plot_pulsed = plt.plot(mySim.tHist, mySim.n2Hist, label = "n2", marker = 's')
plt.xlabel("t")
plt.ylabel("n_i")
plt.legend()
plt.show()

#Do the filename with all the parameters of the simulation   
filename = "PULSER1_HIST_sim_N={0}_r0={1}_r1={2}_r2={3}_u1={4}_u2={5}_SPDP={6}".format(mySim.N, mySim.r0, mySim.r1, mySim.r2, mySim.u1, mySim.u2, simsPerDataPoint)

#Save the data to a file
#SimTools.SaveXYToFile(filename, dataPointsX, dataPointsY)

        
