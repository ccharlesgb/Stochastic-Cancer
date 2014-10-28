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
mySim.timeLimit = 1000.0
mySim.u1 = 0.1
mySim.u2 = 0.1
mySim.r1 = 0.75
mySim.r2 = 0.5
mySim.populationHistory = 1

dataPointCount = 10

simsPerDataPoint = 10

pulseAmp = mySim.r1
pulseWidth = 5.0
maxPulseWidth = 10.0

pulseWavelength = 10.0

dataPulseTime = []
dataPulsePulse = []

FPTHist = []
BoundaryPos0 = []
BoundaryPos1 = []
BoundaryPos2 = []

ReactivePos0 = []
ReactivePos1 = []
ReactivePos2 = []

#Pre Simulate callback (called every frame before a timestep)
def pulse_r1(sim):
    global curPoint
    global dataPointCount
    global pulseWidth

    sim.r1 = SimTools.PulseWave(sim.curTime, pulseAmp, pulseWidth, pulseWavelength) + pulseAmp
    
    dataPulseTime.append(sim.curTime)
    dataPulsePulse.append(sim.r1)  
    
mySim.preSim = pulse_r1 #IMPORTANT assign the callback (called in the class sim loop)

lastFPData = 0.0

def raise_r1(sim):
    global BoundaryPos
    global lastFPData
    r_grad = 0.002
    sim.r1 = r_grad * sim.curTime
    
    dataPulseTime.append(sim.curTime)
    dataPulsePulse.append(sim.r1)
    
    BP = BoundaryFP(sim.r1,sim.r2,sim.u2)
    RP = ReactiveFP(sim.r0,sim.r1,sim.r2, sim.u1,sim.u2)
    if lastFPData + 10.0 < sim.curTime:
        FPTHist.append(sim.curTime)
        BoundaryPos0.append(BP[0] * sim.N)
        BoundaryPos1.append(BP[1] * sim.N)
        BoundaryPos2.append(BP[2] * sim.N)
        
        ReactivePos0.append(RP[0] * sim.N)
        ReactivePos1.append(RP[1] * sim.N)
        ReactivePos2.append(RP[2] * sim.N)
        lastFPData = sim.curTime
    

mySim.preSim = raise_r1

#Run the simulation
mySim.Simulate()

plt.subplot(2,1,1)
plt.plot(dataPulseTime, dataPulsePulse)
plt.xlabel("time")
plt.ylabel("r1")

#Create graph of data
plt.subplot(2, 1, 2)
plot_pulsed = plt.plot(mySim.tHist, mySim.n0Hist, label = "n0", linewidth = 0.25)
plot_pulsed = plt.plot(mySim.tHist, mySim.n1Hist, label = "n1", linewidth = 0.25)
plot_pulsed = plt.plot(mySim.tHist, mySim.n2Hist, label = "n2", linewidth = 0.25)

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

        
