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
mySim.timeLimit = 100.0
mySim.u1 = 0.1
mySim.u2 = 0.1
mySim.r1 = 1.1
mySim.r2 = 1.2
mySim.populationHistory = 0

mySim.in0 = 10
mySim.in1 = 45
mySim.in2 = 45

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
    sim.r1 = SimTools.PulseWave(sim.curTime, pulseOn - pulseOff, pulseWidth, pulseWavelength) + pulseOff
    
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

simsPerDataPoint = 10

#Initialize the array with default values
dataPointsX = []
dataPointsY = []

maxArea = 1.0

wavMin = 0.1
wavMax = 1.0

for curPoint in range(0, dataPointCount):
    startTime = time.clock() #Algorithm benchmarking
    
    fixationCounts = 0
    
    pulseWavelength = float(curPoint)/(dataPointCount-1) * (wavMax - wavMin) + wavMin
    pulseWidth = pulseWavelength / 2.0    
    
    #pulseOn = float(curPoint)/(dataPointCount-1) * (drugOnMax - drugOnMin) + drugOnMin
    pulseOn = maxArea / pulseWidth
    
    print("Pulse on: {0}".format(pulseOn))
    print("Current Data Point = {0}/{1} ({2}%)".format(curPoint + 1, dataPointCount, 100.0 * float(curPoint+1.0)/dataPointCount))
    #Perform many simulations to get an accurate probability of the fixation probability
    for sim in range(0, simsPerDataPoint):
        mySim.Simulate()
            
        if mySim.n2 == mySim.N: #The simulation ended with fixation
            fixationCounts += 1
    #Once the loop is done get the fraction of fixations for this r1
    dataPointsX.append(pulseWavelength)
    dataPointsY.append(float(fixationCounts) / float(simsPerDataPoint))

    print("Complete (Took {:.{s}f} seconds)".format(time.clock() - startTime, s=2))

plot_pulsed = plt.plot(dataPointsX, dataPointsY, label = "Fixation", linewidth = 1.0)
plt.xlabel("Pulse Wavelength")
plt.ylabel("Fixation")
plt.legend()
plt.show()

#Do the filename with all the parameters of the simulation   
filename = "PULSER1_HIST_sim_N={0}_r0={1}_r1={2}_r2={3}_u1={4}_u2={5}_SPDP={6}".format(mySim.N, mySim.r0, mySim.r1, mySim.r2, mySim.u1, mySim.u2, simsPerDataPoint)

#Save the data to a file
#SimTools.SaveXYToFile(filename, dataPointsX, dataPointsY)

        
