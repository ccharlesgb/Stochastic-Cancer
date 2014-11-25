# -*- coding: utf-8 -*-
"""
Created on Sat Nov  8 10:28:27 2014

@author: Jonny
"""
import math
import TwoSpecies
import SimTools
import matplotlib.pyplot as plt
import multipulse
'''
mySim = SimTools.Gillespie(100)
mySim.r0 = 1.0
mySim.r1 = 1.0
mySim.r2 = 1.0
mySim.u1 = 0.1
mySim.u2 = 0.01

strength = 0.6 #drug strength
angle = math.pi/4.0 #angle of pulsing in r1/r2 space
init_r1 = 1.3
init_r2 = 1.3
r1_offset_frac = 0.0
r2_offset_frac = 0.25
freq=1.0/10.0
r1_width_frac=0.5
r2_width_frac=0.5
'''
def multiple_pulse(sim, freq, drug_strength,angle, init_r1, init_r2, r1_offset=0.0, r2_offset=0.0, r1_width=0.5, r2_width=0.5, ret = None):
    #global freq
    global dataPointCount
    global pulseWidth
    global curTime
    global treat_time

    print("angle ARRIVED at multiple pulse".format(angle))
    #print("freq ARRIVED at multiple pulse".format(freq))
    r1_pulseOn = drug_strength*math.cos(angle) + init_r1   #set the maximum and minimum values of the pulses r1
    r1_pulseOff = init_r1    
    
    r2_pulseOn=drug_strength*math.sin(angle) + init_r2
    r2_pulseOff=init_r2    
    
    frequency = freq
    
    sim.r1 = TwoSpecies.PulseWave(sim.curTime, r1_pulseOn - r1_pulseOff, frequency, offset=r1_offset, width_frac=r1_width) + r1_pulseOff #set the simulations r1 to the pulse wave
    sim.r2 = TwoSpecies.PulseWave(sim.curTime, r2_pulseOn - r2_pulseOff, frequency, offset=r2_offset, width_frac=r2_width) + r2_pulseOff
    
    if ret !=None: #optional parameter to allow east plotting of the pulse function
        return (TwoSpecies.PulseWave(sim.curTime, r1_pulseOn - r1_pulseOff, frequency, offset=r1_offset, width_frac=r1_width) + r1_pulseOff , TwoSpecies.PulseWave(sim.curTime, r2_pulseOn - r2_pulseOff, frequency, offset=r2_offset, width_frac=r2_width) + r2_pulseOff)


mySim = SimTools.Gillespie(10)

mySim.in0 = 10
mySim.timeLimit = 4000.0
mySim.u1 = 0.1
mySim.u2 = 0.1
mySim.r1 = 1.1
mySim.r2 = 1.3
mySim.populationHistory = 0

mySim.in0 = 10
mySim.in1 = 0
mySim.in2 = 0

#pick initial r1, to find the starting point in r1/r2 space
r1_origin=1.25
r2_origin=(1.0-mySim.u2)*r1_origin
mut2= math.pow((1.0-mySim.u2),2) #define a term that appears alot in finding the initial r1 and r2

#define frequency for both waves
frequency = 0.05

#set up multiple pulse parameters
myPulse= multipulse.Multipulse()

myPulse.drug_strength = 0.6
myPulse.angle = math.pi/2.0

myPulse.init_r1 = -math.pow( ( mut2 / (1.0+mut2) ), 0.5) * myPulse.drug_strength/5.0 + r1_origin
myPulse.init_r2 = math.pow( (1.0/(1.0+mut2) ) ,0.5)* myPulse.drug_strength/5.0 + r2_origin        
        
#define r1 pulse paramters
myPulse.freq_r1 = frequency
myPulse.r1_width = 0.1
myPulse.r1_offset = 0.0
        
#define r2 pulse parameters       
myPulse.freq_r2 = frequency
myPulse.r2_width = 0.1
myPulse.r2_offset = 0.0


mySim.pulseParam=myPulse

mySim.preSim=multipulse.multiple_pulse


xmax = 400
r1_values=[]
r2_values=[]
curTime = range(0,int(xmax))  
for j in curTime:
    curPoint=j
    mySim.curTime=curTime[j]
    r1_values.append(mySim.r1)
    r2_values.append(mySim.r1)




   
plt.subplot(2, 1, 1)
plt.plot(curTime,r1_values)
plt.xlabel("Time")
plt.ylabel("r1")

plt.subplot(2, 1, 2)
plt.plot(curTime,r2_values)
plt.xlabel("Time")
plt.ylabel("r2")

'''
plt.subplot(3, 1, 3)
plt.plot(curTime,total_pulse)
plt.xlabel("Time")
plt.ylabel("Total")
plt.show()
'''