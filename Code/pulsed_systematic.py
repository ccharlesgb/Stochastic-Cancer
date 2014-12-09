# -*- coding: utf-8 -*-
"""
Created on Thu Nov 27 17:44:42 2014

@author: Jonny
"""

import systematic_transitions
import SimTools
import numpy as np
import math
import matplotlib.pyplot as plt
import MatTools
import TwoSpecies

#general simulation parameters
mySim = SimTools.Gillespie(10)

mySim.in0 = 10
mySim.timeLimit = 16000.0
mySim.u1 = 0.1
mySim.u2 = 0.1
mySim.r0 = 1.0
mySim.r1 = 1.25
mySim.r2 = 1.125

mySim.ResetSim()

#assign some pulse parameters
max_freq = 0.5
min_freq = 0.1

#create a pulsed function for type-1 fitness
def pulse_r1_num(sim):
    global curPoint
    global DPC
    global pulseWidth
    global curTime
    global multiplier
    
    pulseOn = 1.3    #set the maximum and minimum values of the pulses r1
    pulseOff = 0.7    
    
    #pulsePeriod = multiplier*((float(dataPointCount) ) / (float(curPoint)+1.0)) * maxPulsePeriod #how the time period varies with time - time period is reciprical to the frequency
    frequency = float(curPoint)/(float(DPC)-1.0)*(max_freq - min_freq) + min_freq

    numer.sim.r1 = TwoSpecies.PulseWave(numer.sim.curTime, pulseOn - pulseOff, frequency) + pulseOff #set the simulations r1 to the pulse wave

def pulse_r1_sim(sim):
    global curPoint
    global DPC
    global pulseWidth
    global curTime
    global multiplier
    
    pulseOn = 1.3    #set the maximum and minimum values of the pulses r1
    pulseOff = 0.7    
    
    #pulsePeriod = multiplier*((float(dataPointCount) ) / (float(curPoint)+1.0)) * maxPulsePeriod #how the time period varies with time - time period is reciprical to the frequency
    frequency = float(curPoint)/(float(DPC)-1.0)*(max_freq - min_freq) + min_freq

    sim.r1 = TwoSpecies.PulseWave(sim.curTime, pulseOn - pulseOff, frequency) + pulseOff #set the simulations r1 to the pulse wave



#assign simulation to callback to r1 function
mySim.preSim = pulse_r1_sim




DPC = 15
SPD = 5

fixTime_sim = []

min_r1 = 0.6
max_r1 = 1.4
r1 = []
errors = []

frequency =[]

mySim.ResetSim()
for curPoint in range(0,DPC): #simulate some shit
    frequency.append(float(curPoint)/(float(DPC)-1.0)*(max_freq - min_freq) + min_freq) #create array to plot frequency    
    print("Data Point {0}/{1} - {2}%".format(curPoint + 1, DPC, 100*(curPoint + 1)/DPC))    
    fixTimeTerm = 0.0    
    for i in range(0,SPD):
        mySim.Simulate()
        if mySim.Fixated():
            fixTimeTerm += mySim.curTime
        else:
            fixTimeTerm += mySim.timeLimit
    fixTime_sim.append(fixTimeTerm/SPD)
    errors.append(fixTime_sim[curPoint]*math.pow(SPD,-0.5))


fixTime_num = []
frequency = []
#
numer = systematic_transitions.systematic(mySim)
 
#assign numerical to callback to r1 function
mySim.preSim = pulse_r1_num
 
 
for curPoint in range(0,DPC):
    frequency.append(float(curPoint)/(float(DPC)-1.0)*(max_freq - min_freq) + min_freq) #create array to plot frequency
    mySim.ResetSim()    
    result = numer.Get_fix_time()
    fixTime_num.append(result)


#MatTools.SaveXYData("Jnumerical_fixation_time_varying_r1_02_to_30", r1, fixTime_num, xLabel = "r1", yLabel = "Fixation Time", sim = mySim)
#MatTools.SaveXYData("Jsimulation_fixation_time_varying_r1_02_30", r1, fixTime_num, yError = errors ,xLabel = "r1", yLabel = "Fixation Time", sim = mySim)

plt.errorbar(frequency,fixTime_sim,errors, label = "Simulation")
plt.plot(frequency,fixTime_num,label = "Numerical", marker = "o")
plt.xlabel("Frequency of r1 pulse")
plt.ylabel("Fixation Time")
plt.legend()
plt.show()







