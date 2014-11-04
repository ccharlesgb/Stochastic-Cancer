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
#import anfixtime
import numpy as np
import copy

#Initialize the Gillespie simulator with 10 cells
N=10
mySim = TwoSpecies.Gillespie(N)
mySim.timeLimit = 3000
dataPointCount = N
#mySim.j = int(N/2)
mySim.r0=1.0
mySim.r1=1.0


pulseOff = 0.5
pulseOn = 1.5
pulseWidth = 1.5
maxPulseWidth = 1.0

maxPulseWavelength = 10.0

#Pre Simulate callback (called every frame before a timestep)

#set up pulsed function for r1
def pulse_r1(sim):
    global curPoint
    global dataPointCount
    global pulseWidth
    global curTime    
    
    pulseWavelength = (float(curPoint) / (dataPointCount - 1)) * maxPulseWavelength
    sim.r1 = SimTools.PulseWave(sim.curTime, pulseOn - pulseOff, pulseWidth, pulseWavelength) + pulseOff


mySim.preSim = pulse_r1 #IMPORTANT assign the callback (called in the class sim loop)
'''
#define average pulse for r1
def avgpulse_r1(sim):
    global curPoint
    global dataPointCount
    global pulseWidth
    pulseWidth = (float(curPoint) / (dataPointCount - 1)) * maxPulseWidth
    sim.r1 = (pulseOn - pulseOff) * (pulseWidth / pulseWavelength) + pulseOff
   ''' 

#compare the variance of the distributions with the pulse time
def GetStandardDev(distribution):
    dist2=[]    
    for i in range(0,len(distribution)):    
        dist2.append(math.pow(distribution[i],2))
    
    term1=sum(dist2)/len(dist2)
    term2=math.pow((sum(distribution)/len(distribution)),2)
    variance = term1-term2
    standard_dev=math.pow(variance,0.5)
    return standard_dev
    


#Sweep the parameter r1 from 0.2 to 3.0 and run many simulations per data point
#Gets an idea on how likely cancer fixation is to occur for this parameter
simsPerDataPoint = 10000

#Initialize the array with default values
dataPointsX = []
dataPointsY = []
for i in range(0, dataPointCount):
    dataPointsX.insert(i,0.0)
    dataPointsY.insert(i,0.0)


#initialise array to count indiviual fixation times to get distribution
    indiv_distrib = []
    total_distrib =[]
#preallocate distribution array for individual simulation
for i in range(0, simsPerDataPoint):
    indiv_distrib.insert(i,0.0)    
#preallocate distribution array for all simulations
for i in range(0,dataPointCount):    
    total_distrib.insert(i,0.0)

'''
#for pltting r1 pulse
pulse_time=[]
pulse_mag=[]
#preallocate array to plot the r1 pulse function
for j in range (0,dataPointCount):
    pulse_time.append(0.0) 
    pulse_mag.append(0.0)
    for i in range(0,mySim.timeLimit):
        pulse_time[j].append(0.0)
        pulse_mag[j].append(0.0)
'''

for curPoint in range(0, dataPointCount):
    startTime = time.clock() #Algorithm benchmarking
    
    fixationTime = 0
    mySim.ij = int(float(curPoint)/(dataPointCount-1) * mySim.N)
    print("Current Data Point = {0}/{1} ({2}%)".format(curPoint + 1, dataPointCount, 100.0 * float(curPoint+1.0)/dataPointCount))
    #Perform many simulations to get an accurate probability of the fixation probability
    for i in range(0, simsPerDataPoint):
        mySim.Simulate()
            
        if mySim.j == mySim.N or mySim.j==0: #The simulation ended with fixation
            fixationTime += mySim.curTime
            indiv_distrib[i] = copy.copy(mySim.curTime) #add distribution of fix times to array for each iteration
            
        
    total_distrib[curPoint]= copy.copy(indiv_distrib) #add each data points fix time distribution to an array containing all distribs
    
    #Once the loop is done get the fraction of fixations for this r1
    dataPointsX[curPoint] = mySim.ij
    dataPointsY[curPoint] = float(fixationTime) / float(simsPerDataPoint)

    print("Complete (Took {:.{s}f} seconds)".format(time.clock() - startTime, s=2))

#plot distribution of fixation times for different data points
for i in range(0, dataPointCount):
    hist, bins = np.histogram(total_distrib[i], bins=50)
    width = (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    plt.bar(center, hist, align='center', width=width)
    plt.title("Data Point:{0} Average Value:{1}".format(i+1,sum(total_distrib[i])/len(total_distrib[i])))    
    plt.xlabel("Fixation Time: ")
    plt.ylabel("Counts")
    plt.ylim([0,7500])
    plt.xlim([0,5000])
    fname="FixTimeDistribFiguresDataPoint_{0}_Frequency={1}MaxPulse={2}MinPulse={3}PulseWidth={4}.png".format(i, pulseWavelength, pulseOn, pulseOff, (float(i) / (dataPointCount - 1)) * maxPulseWidth)
    plt.savefig(fname)
    plt.show()

standard_dev=[]
#get the variance against thepulse time
for i in range(0,dataPointCount):
    standard_dev.append(GetStandardDev(total_distrib[i]))

plot_standard_dev= plt.plot(dataPointsX, standard_dev, label="Standard Dev")
plt.title("Standard Deviation vs r1 Pulse Time")
plt.xlabel("r1 Pulse Time")
plt.ylabel("Standard Dev")
plt.show()


#Dump data to console
for i in range(0, dataPointCount):
    print("radTime: {0} Fixation: {1}".format(dataPointsX[i],dataPointsY[i]))

#Create graph of data
plot_pulsed = plt.errorbar(dataPointsX, dataPointsY, math.pow(simsPerDataPoint,-0.5), label = "Pulsed")
plt.xlabel("r1 Pulse Time: ")
plt.ylabel("Fixation Time")
plt.show()

'''
#initial crude way to find local maxima, and hence difference between peaks
def FindLocalMaxima(array):
    maxima=[]    
    for i in range(0, len(array)):
        if array[i]>array[i-1] and array[i]>array[i+1]:
            maxima.append[i]
    if len(maxima) == 2:
        return maxima[1]-maxima[0]
    if len(maxima) !=2:
        print("There is more than two local maximum. Panic.")
        return 0

difference=[]
for i in range(0, dataPointCount):
    difference.append(FindLocalMaxima(total_distrib[i]))
    
#plot peak difference against 
datapointsX=[]
datapointsY=[]
peaks=[]
for i in range(0,dataPointCount):
    datapointsX.append(0.0)
    datapointsY.append(0.0)
    peaks.append(0.0)
    
for i in range(0,dataPointCount):
    peaks[i]=FindLocalMaxima(total_distrib[i])



plot_pulsed = plt.plot(dataPointsX, datapointsX, label = "Pulsed")
plt.xlabel("r1 Pulse Time: ")
plt.ylabel("peak diff")
plt.show()
'''
#Do the filename with all the parameters of the simulation   
filename = "TwoSpeciesPulsed_sim_N={0}_r0={1}_r1={2}_SPDP={3}".format(mySim.N, mySim.r0, mySim.r1, simsPerDataPoint)

#Save the data to a file
#SimTools.SaveXYToFile(filename, dataPointsX, dataPointsY)

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

plot_averaged = plt.errorbar(dataPointsX, dataPointsY, math.pow(simsPerDataPoint,0.5), label = "Averaged")

plt.legend()

plt.show()


#theoretical non-conditional fixation time (for averaged p)
dataPointsX = []
dataPointsY = []
for i in range(0,mySim.N):
    dataPointsX.append(i,i)
    dataPointsY.append(i,anfixtime.GetFixTimeJ(sim,i))

plt.plot(dataPointsX, dataPointsY, label = "Theoretical")'''