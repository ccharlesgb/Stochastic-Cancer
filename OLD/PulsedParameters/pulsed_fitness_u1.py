import matplotlib.pyplot as plt
import SimTools
import time
import math

import Vi

print(Vi.GetV_i(1, 1.0,1.0,1.0, 0.1, 0.01, 10))

#Initialize the Gillespie simulator with 10 cells
mySim = SimTools.Gillespie(10)
mySim.timeLimit = 100.0
mySim.u1 = 0.1
mySim.u2 = 0.01
mySim.r1 = 1.0
mySim.r2 = 1.0
dataPointCount = 10

pulseOff = 0.01
pulseOn = 0.1
pulseWidth = 1.0
maxPulseWidth = 10.0

pulseWavelength = 10.0

#Pre Simulate callback (called every frame before a timestep)
def pulse_u1(sim):
    global curPoint
    global dataPointCount
    global pulseWidth

    pulseWidth = (float(curPoint) / (dataPointCount - 1)) * maxPulseWidth
    sim.u1 = SimTools.PulseWave(sim.curTime, pulseOn - pulseOff, pulseWidth, pulseWavelength) + pulseOff

mySim.preSim = pulse_u1 #IMPORTANT assign the callback (called in the class sim loop)

def avgpulse_u1(sim):
    global curPoint
    global dataPointCount
    global pulseWidth
    pulseWidth = (float(curPoint) / (dataPointCount - 1)) * maxPulseWidth
    sim.u1 = (pulseOn - pulseOff) * (pulseWidth / pulseWavelength) + pulseOff
    

#Sweep the parameter r1 from 0.2 to 3.0 and run many simulations per data point
#Gets an idea on how likely cancer fixation is to occur for this parameter
simsPerDataPoint = 20000
fixError = 1/math.sqrt(float(simsPerDataPoint))

#Initialize the array with default values
dataPointsX = []
dataPointsY = []
for i in range(0, dataPointCount):
    dataPointsX.insert(i,0.0)
    dataPointsY.insert(i,0.0)

for curPoint in range(0, dataPointCount):
    startTime = time.clock() #Algorithm benchmarking
    
    fixationCounts = 0
    
    print("Current Data Point = {0}/{1} ({2}%)".format(curPoint + 1, dataPointCount, 100.0 * float(curPoint+1.0)/dataPointCount))
    #Perform many simulations to get an accurate probability of the fixation probability
    for sim in range(0, simsPerDataPoint):
        mySim.Simulate()
            
        if mySim.n2 == mySim.N: #The simulation ended with fixation
            fixationCounts += 1
    #Once the loop is done get the fraction of fixations for this r1
    dataPointsX[curPoint] = pulseWidth
    dataPointsY[curPoint] = float(fixationCounts) / float(simsPerDataPoint)

    print("Complete (Took {:.{s}f} seconds)".format(time.clock() - startTime, s=2))

#Dump data to console
for i in range(0, dataPointCount):
    print("radTime: {0} Fixation: {1}".format(dataPointsX[i],dataPointsY[i]))

#Create graph of data
plot_pulsed = plt.errorbar(dataPointsX, dataPointsY, yerr = fixError, label = "Pulsed")
plt.xlabel("u1 Pulse Time: ")
plt.ylabel("Type 2 Fixation Prob")
plt.show()

#Do the filename with all the parameters of the simulation   
filename = "PULSEU1_sim_N={0}_r0={1}_r1={2}_r2={3}_u1={4}_u2={5}_SPDP={6}".format(mySim.N, mySim.r0, mySim.r1, mySim.r2, mySim.u1, mySim.u2, simsPerDataPoint)

#Save the data to a file
SimTools.SaveXYToFile(filename, dataPointsX, dataPointsY)

#Do the average of the pulse and see how it compares
#Initialize the array with default values

mySim.preSim = avgpulse_u1

dataPointsX = []
dataPointsY = []
for i in range(0, dataPointCount):
    dataPointsX.insert(i,0.0)
    dataPointsY.insert(i,0.0)

for curPoint in range(0, dataPointCount):
    startTime = time.clock() #Algorithm benchmarking
    
    fixationCounts = 0
    
    print("Current Data Point = {0}/{1} ({2}%)".format(curPoint + 1, dataPointCount, 100.0 * float(curPoint+1.0)/dataPointCount))
    #Perform many simulations to get an accurate probability of the fixation probability
    for sim in range(0, simsPerDataPoint):
        mySim.Simulate()
            
        if mySim.n2 == mySim.N: #The simulation ended with fixation
            fixationCounts += 1
    #Once the loop is done get the fraction of fixations for this r1
    dataPointsX[curPoint] = pulseWidth
    dataPointsY[curPoint] = float(fixationCounts) / float(simsPerDataPoint)

    print("Complete (Took {:.{s}f} seconds)".format(time.clock() - startTime, s=2))

#Dump data to console
for i in range(0, dataPointCount):
    print("radTime: {0} Fixation: {1}".format(dataPointsX[i],dataPointsY[i]))

plot_averaged = plt.errorbar(dataPointsX, dataPointsY, yerr = fixError, label = "Averaged")

plt.legend()

plt.show()


        
