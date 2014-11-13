# -*- coding: utf-8 -*-
"""
Created on Sat Nov  8 10:28:27 2014

@author: Jonny
"""
import math
import TwoSpecies
import SimTools
import matplotlib.pyplot as plt
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
'''
    if plot != None:
       curTime = range(0,int(sim.timeLimit))
       r1_pulse=[]
       r2_pulse=[]
       total_pulse=[]
       for j in curTime:
           curPoint=j
           sim.curTime=curTime[j]
           r1_pulse_term , r2_pulse_term = multiple_pulse(sim, drug_strength, angle, init_r1, init_r2, r1_offset=r1_offset, r2_offset=r2_offset, r1_width=r1_width , r2_width=r2_width, ret=1)
           r1_pulse.append(r1_pulse_term)
           r2_pulse.append(r2_pulse_term)
           total_pulse.append(math.pow((math.pow(r1_pulse_term,2)+math.pow(r2_pulse_term,2)),0.5))
       plt.subplot(3, 1, 1)
       plt.plot(curTime,r1_pulse)
       plt.xlabel("Time")
       plt.ylabel("r1")
       plt.title
       plt.subplot(3, 1, 2)
       plt.plot(curTime,r2_pulse)
       plt.xlabel("Time")
       plt.ylabel("r2")
       plt.subplot(3, 1, 3)
       plt.plot(curTime,total_pulse)
       plt.xlabel("Time")
       plt.ylabel("Total")
       plt.title(plot)
       plt.show()
       
       
#mySim.preSim=multiple_pulse #assign the multiple pulse function to the simulation
'''

def plot_multi_pulse(sim, drug_strength, angle, init_r1, init_r2, r1_offset=0.0, r2_offset=0.0, r1_width=0.5, r2_width=0.5, Title=None):
    curTime = range(0,40)
    r1_pulse=[]
    r2_pulse=[]
    total_pulse=[]    
    
    for j in curTime:
       curPoint=j
       sim.curTime=curTime[j]
       r1_pulse_term , r2_pulse_term = multiple_pulse(sim, drug_strength, angle, init_r1, init_r2, r1_offset=r1_offset, r2_offset=r2_offset, r1_width=r1_width, r2_width=r2_width, ret=1)
       r1_pulse.append(r1_pulse_term)
       r2_pulse.append(r2_pulse_term)
       total_pulse.append(math.pow((math.pow(r1_pulse_term,2)+math.pow(r2_pulse_term,2)),0.5))
       
    plt.subplot(3, 1, 1)
    plt.plot(curTime,r1_pulse)
    plt.xlabel("Time")
    plt.ylabel("r1")
    plt.title
    plt.subplot(3, 1, 2)
    plt.plot(curTime,r2_pulse)
    plt.xlabel("Time")
    plt.ylabel("r2")
    
    plt.subplot(3, 1, 3)
    plt.plot(curTime,total_pulse)
    plt.xlabel("Time")
    plt.ylabel("Total")
    plt.show()
    plt.title(Title)