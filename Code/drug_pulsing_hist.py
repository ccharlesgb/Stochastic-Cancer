# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import SimTools
import time
import math
import numpy

def FixedPointExists(r1, r2, u2):
    condition = r2 / (1.0 - u2)
    return r1 > condition and r2 < condition
    
def BoundaryFP(r1, r2, u2):
    exists = (r1 > r2/(1.0-u2)) and (r2 < r2/(1.0-u2))
    if exists == 0:
        fp = []
        fp.append(numpy.nan)
        fp.append(numpy.nan)
        fp.append(numpy.nan)
        return fp
    fp = []
    fp.append(numpy.nan) #Append n0
    fp.append(1.0 - u2*r1/(r1-r2))
    fp.append(u2*r1/(r1-r2))
    return fp
    
def ReactiveFP(r0, r1, r2, u1, u2):
    denom = u2*r1*(r0-r2) + (r0-r1)*((1.0-u1)*r0 - r2)
    exists = ((1.0-u1)*r0 > (1.0-u2)*r1) and ((1.0-u1)*r0 > r2)
    if exists == 0:
        fp = []
        fp.append(numpy.nan)
        fp.append(numpy.nan)
        fp.append(numpy.nan)
        return fp
    numer1 = ((1.0-u1)*r0 + u2*r1 - r2)*u1*r0
    fp = []
    fp.append(1.0 - numer1/denom)
    fp.append((((1.0-u1)*r0-r2)*u1*r0)/denom)
    fp.append(u1*u2*r0*r1/denom)
    return fp
    
def FixedPointStable(r0, r1, u2):
    return (1.0 - u2)*r1 > (1.0 - u2) * r0
def FixedPointSaddle(r0,r1, u1, u2):
    exists = (1.0 - u2)*r1 < (1.0 - u1)*r0
    if exists == 0:
        return -1.0
    return 

#Initialize the Gillespie simulator with 10 cells
mySim = SimTools.Gillespie(100)
mySim.timeLimit = 500.0
mySim.u1 = 0.01
mySim.u2 = 0.01
mySim.r1 = 1.1
mySim.r2 = 1.2
mySim.populationHistory = 1

mySim.in0 = 100
mySim.in1 = 0
mySim.in2 = 0

dataPointCount = 10

simsPerDataPoint = 10

fpRecordInterval = 1.0
FPTHist = []
BoundaryPos0 = []
BoundaryPos1 = []
BoundaryPos2 = []

ReactivePos0 = []
ReactivePos1 = []
ReactivePos2 = []

#Pulse Parameters
dataPulseTime = []
dataPulsePulse = []

n0Threshold = 0.2
treatmentOn = 0

pulseOff = mySim.r1
pulseOn = 0.8
pulseWidth = 5.0
maxPulseWidth = 10.0

pulseWavelength = 10.0

#Pre Simulate callback (called every frame before a timestep)
def pulse_r2(sim):
    global treatmentOn
    if float(sim.n0) / sim.N > n0Threshold and treatmentOn == 0:
        sim.r1 = pulseOff
    else:
        treatmentOn = 1
        sim.r1 = SimTools.PulseWave(sim.curTime, pulseOn - pulseOff, pulseWidth, pulseWavelength) + pulseOff
        sim.r2 = sim.r1 + 0.1
    
    dataPulseTime.append(sim.curTime)
    dataPulsePulse.append(sim.r2)    
    
    RecordFP(sim)

mySim.preSim = pulse_r2 #IMPORTANT assign the callback (called in the class sim loop)

lastFPData = 0.0

def RecordFP(sim):
    global lastFPData
    BP = BoundaryFP(sim.r1,sim.r2,sim.u2)
    RP = ReactiveFP(sim.r0,sim.r1,sim.r2, sim.u1,sim.u2)
    if lastFPData + fpRecordInterval < sim.curTime:
        FPTHist.append(sim.curTime)
        BoundaryPos0.append(BP[0] * sim.N)
        BoundaryPos1.append(BP[1] * sim.N)
        BoundaryPos2.append(BP[2] * sim.N)
        
        ReactivePos0.append(RP[0] * sim.N)
        ReactivePos1.append(RP[1] * sim.N)
        ReactivePos2.append(RP[2] * sim.N)
        lastFPData = sim.curTime

#Run the simulation
mySim.Simulate()

plt.subplot(2,1,1)
plt.plot(dataPulseTime, dataPulsePulse)
plt.xlabel("time")
plt.ylabel("r2")

#Create graph of data
plt.subplot(2, 1, 2)
plot_pulsed = plt.plot(mySim.tHist, mySim.n0Hist, label = "n0", linewidth = 0.5)
plot_pulsed = plt.plot(mySim.tHist, mySim.n1Hist, label = "n1", linewidth = 0.5)
plot_pulsed = plt.plot(mySim.tHist, mySim.n2Hist, label = "n2", linewidth = 0.5)

plot_pulsed = plt.plot(FPTHist, BoundaryPos0, label = "Boundary0", marker = 'o', linestyle = '')
plot_pulsed = plt.plot(FPTHist, BoundaryPos1, label = "Boundary1", marker = '^', linestyle = '')
plot_pulsed = plt.plot(FPTHist, BoundaryPos2, label = "Boundary2", marker = 's', linestyle = '')

plot_pulsed = plt.plot(FPTHist, ReactivePos0, label = "Reactive1", marker = 'o', linestyle = '')
plot_pulsed = plt.plot(FPTHist, ReactivePos1, label = "Reactive2", marker = '^', linestyle = '')
plot_pulsed = plt.plot(FPTHist, ReactivePos2, label = "Reactive3", marker = 's', linestyle = '')
plt.xlabel("t")
plt.ylabel("n_i")
plt.legend()
plt.show()

#Do the filename with all the parameters of the simulation   
filename = "PULSER1_HIST_sim_N={0}_r0={1}_r1={2}_r2={3}_u1={4}_u2={5}_SPDP={6}".format(mySim.N, mySim.r0, mySim.r1, mySim.r2, mySim.u1, mySim.u2, simsPerDataPoint)

#Save the data to a file
#SimTools.SaveXYToFile(filename, dataPointsX, dataPointsY)

        
