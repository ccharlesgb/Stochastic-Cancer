# -*- coding: utf-8 -*-
"""
Created on Mon Nov 24 17:37:39 2014

@author: Jonny
"""

import matplotlib.pyplot as plt
import SimTools
import time
import math
import numpy
import multipulse



#set up simulation parameters
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

#define simulation parameters
DPC = 20
SPD = 25000

min_offset=0.0
max_offset=1.0
fixTime=[]
phase=[]
errors=[]
for curPoint in range(0,DPC):
    offset=(curPoint/(DPC-1.0))*(max_offset-min_offset) + min_offset
    phase.append(offset*2*math.pi)    
    myPulse.r2_offset = offset
    fixTime_term = 0.0
    print("Current Data Point = {0}/{1} ({2}%)".format(curPoint + 1.0, DPC, 100.0 * float(curPoint+1.0)/DPC))
    for i in range (0, SPD):
        mySim.Simulate()
        if mySim.Fixated():
            fixTime_term += mySim.curTime
        else:
            fixTime_term += mySim.timeLimit
            print("WARNING! - system did not fixate")
    average_fix=fixTime_term/SPD
    fixTime.append(average_fix)
    errors.append(average_fix*math.pow(SPD,-0.5))

plt.errorbar(phase, fixTime, errors)
plt.xlabel("Phase Difference")
plt.ylabel("Fixation Time")
plt.title("Phase Difference vs Fixation Time")
fname="Phase_difference_vs_fixtime"
plt.savefig(fname)























