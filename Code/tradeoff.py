# -*- coding: utf-8 -*-
"""
Created on Thu Nov 13 16:05:50 2014

@author: Jonny
"""

import matplotlib.pyplot as plt
import TwoSpecies
import time
import math
#import anfixtime
import numpy as np
import copy
import multipulse

#Initialize the Gillespie simulator with 10 cells
N=10
mySim = TwoSpecies.Gillespie(N)
mySim.timeLimit = 20000
dataPointCount = 2
mySim.r0=1.0
mySim.r1=1.0
mySim.u1=0.1

simsPerDataPoint = 10000
dataPointCount = 50

min_freq = 0.01
max_freq = 0.5


#
def pulse_r1(sim):
    global curPoint
    global dataPointCount
    global pulseWidth
    global curTime
    global multiplier
    
    pulseOn = 1.0    #set the maximum and minimum values of the pulses r1
    pulseOff = 0.5    
    
    #pulsePeriod = multiplier*((float(dataPointCount) ) / (float(curPoint)+1.0)) * maxPulsePeriod #how the time period varies with time - time period is reciprical to the frequency
    frequency = float(curPoint)/(float(dataPointCount)-1.0)*(max_freq-min_freq) + min_freq
    #pulseWidth=pulsePeriod/2.0  #this is to keep the average value constatnt
    
    sim.r1 = TwoSpecies.PulseWave(sim.curTime, pulseOn - pulseOff, frequency, offset = 0.5) + pulseOff #set the simulations r1 to the pulse wave

def pulse_u1(sim):
    global curPoint
    global dataPointCount
    global pulseWidth
    global curTime
    global multiplier
    
    pulseOn = 0.4    #set the maximum and minimum values of the pulses u1
    pulseOff = 0.1
    
    #pulsePeriod = multiplier*((float(dataPointCount) ) / (float(curPoint)+1.0)) * maxPulsePeriod #how the time period varies with time - time period is reciprical to the frequency
    frequency = float(curPoint)/(float(dataPointCount)-1.0)*(max_freq-min_freq) + min_freq
    #pulseWidth=pulsePeriod/2.0  #this is to keep the average value constatnt
    
    sim.u1 = TwoSpecies.PulseWave(sim.curTime, pulseOn - pulseOff, frequency, offset = 0.0) + pulseOff #set the simulations r1 to the pulse wave
            

def dummy_callback(sim):
    pulse_r1(sim)    
    pulse_u1(sim)

#assign the callback
mySim.preSim=dummy_callback
fixTime = []
frequency=[]
for curPoint in range (0,dataPointCount):
    frequency.append((curPoint/(dataPointCount + 1.0))*(max_freq-min_freq) + min_freq)

    fixTimeTerm=0.0
    for j in range (0,simsPerDataPoint):
        mySim.Simulate()
        if mySim.Fixated():
            fixTimeTerm+=mySim.curTime
        else:
            fixTimeTerm += mySim.timeLimit
            print("Warning! System did not fixate")
    
    fixTime.append(fixTimeTerm/simsPerDataPoint)
    
plt.plot(frequency, fixTime)
plt.xlabel("Frequency")
plt.ylabel("Fixation Time")
plt.show()

        
        
        
        
        
        
        
        
        
        
        
        
        

