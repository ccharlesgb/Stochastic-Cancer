# -*- coding: utf-8 -*-
"""
Created on Thu Oct 23 11:41:06 2014

@author: Jonny
"""

import matplotlib.pyplot as plt
import TwoSpecies
import time
import math
#import anfixtime
import numpy as np
import copy

#Initialize the Gillespie simulator with 10 cells
N=10
mySim = TwoSpecies.Gillespie(N)
mySim.timeLimit = 3000
dataPointCount = 50
mySim.r0=1.0
mySim.r1=1.0
mySim.u1=0.1

simsPerDataPoint = 5000
dataPointCount = 25

#average_r1=1.0
#pulseWidth = 2.0
#maxPulseWidth = 1.0
#pulseOn=1.5
#pulseOff=0.5

maxPulsePeriod = 2.0
multiplier=10.0
maxFrequency = 1.0/28.0

#Pre Simulate callback (called every frame before a timestep)

#set up pulsed function for r1
def pulse_r1(sim, ret = None):
    global curPoint
    global dataPointCount
    global pulseWidth
    global curTime
    global multiplier
    
    pulseOn = 1.3    #set the maximum and minimum values of the pulses r1
    pulseOff = 0.7    
    
    #pulsePeriod = multiplier*((float(dataPointCount) ) / (float(curPoint)+1.0)) * maxPulsePeriod #how the time period varies with time - time period is reciprical to the frequency
    frequency = float(curPoint)/(float(dataPointCount)-1.0)*maxFrequency
    #pulseWidth=pulsePeriod/2.0  #this is to keep the average value constatnt
    
    sim.r1 = TwoSpecies.PulseWave(sim.curTime, pulseOn - pulseOff, frequency) + pulseOff #set the simulations r1 to the pulse wave
    
    if ret !=None: #optional parameter to allow east plotting of the pulse function
        return TwoSpecies.PulseWave(sim.curTime, pulseOn - pulseOff, frequency) + pulseOff

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

#defien function to get the standard deviation of a distribution
def GetStandardDev(distribution):
    dist2=[]    #intitalise the squared distribution array
    for i in range(0,len(distribution)):
        dist2.append(math.pow(distribution[i],2)) #square each term in the distribution
    
    term1=sum(dist2)/len(dist2) #find <x^2>
    term2=math.pow((sum(distribution)/len(distribution)),2) #find <x>^2
    variance = term1-term2 #get the variance by subtracting the two terms
    standard_dev=math.pow(variance,0.5) #root the variance to get the standard deviation
    return standard_dev
    


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

#SIMULATION
for curPoint in range(0, dataPointCount):
    startTime = time.clock() #Algorithm benchmarking
    
    fixationTime = 0
    #successful_fixations =  0 #to exclude simulations that fix at j=0 
    #mySim.ij = int(float(curPoint)/(dataPointCount-1) * mySim.N)
    mySim.ij=1
    print("Current Data Point = {0}/{1} ({2}%)".format(curPoint + 1, dataPointCount, 100.0 * float(curPoint+1.0)/dataPointCount))
    
    
    #Perform many simulations to get an accurate probability of the fixation probabilities
    #while successful_fixations < simsPerDataPoint:    
    for i in range(0, simsPerDataPoint):
        mySim.Simulate()
            
        if mySim.j == mySim.N: #or mySim.j==0: #The simulation ended with fixation
            fixationTime += mySim.curTime
            indiv_distrib[i] = copy.copy(mySim.curTime) #add distribution of fix times to array for each iteration
            #successful_fixations += 1 
    
    total_distrib[curPoint]= copy.copy(indiv_distrib) #add each data points fix time distribution to an array containing all distribs
    
    #Once the loop is done get the fraction of fixations for this r1
    dataPointsX[curPoint] = float(curPoint)/(float(dataPointCount)-1.0)*maxFrequency #this is the pulse frequency
    dataPointsY[curPoint] = float(fixationTime) / (float(simsPerDataPoint)) #fixation time
    print("Complete (Took {:.{s}f} seconds)".format(time.clock() - startTime, s=2))

    

#plot distribution of fixation times for different data points
for i in range(0, dataPointCount):
    hist, bins = np.histogram(total_distrib[i], bins=50) #create histogram for each data points distribution with set number of bins
    width = (bins[1] - bins[0]) 
    center = (bins[:-1] + bins[1:]) / 2
    
    plt.subplot(2, 1, 1)    
    plt.bar(center, hist, align='center', width=width)
    plt.title("Data Point:{0} Average Value:{1}".format(i+1,sum(total_distrib[i])/len(total_distrib[i])))    
    plt.xlabel("Fixation Time")
    plt.ylabel("Counts")
    #fname="FixTimeDistribFiguresDataPoint_{0}_Frequency={1}MaxPulse={2}MinPulse={3}PulseWidth={4}.png".format(i, pulseWavelength, pulseOn, pulseOff, (float(i) / (dataPointCount - 1)) * maxPulseWidth)
    #plt.savefig(fname)
    
    # create subplot of the pulse   
    xmin, xmax = plt.xlim()#get the value of the 'maximum' fix time on the simulation  
    y_values=[]
    curTime = range(0,int(xmax))  
    for j in curTime:
        curPoint=i
        mySim.curTime=curTime[j]
        y_values.append(pulse_r1(mySim, 1))
    
    plt.subplot(2, 1, 2)
    plt.plot(curTime,y_values)
    
    
    
    
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
plt.xlabel("r1 Pulse Frequency: ")
plt.ylabel("Fixation Time")
plt.show()


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
'''