# -*- coding: utf-8 -*-
"""
Created on Tue Nov 11 13:26:36 2014

@author: Jonny
"""

import matplotlib.pyplot as plt
import SimTools
import time
import math
import numpy
import multiple_pulse
        


myPulseParam = PulseParam()


#set up simulation parameters
mySim = SimTools.Gillespie(10)

mySim.in0 = 10
mySim.timeLimit = 500.0
mySim.u1 = 0.1
mySim.u2 = 0.1
mySim.r1 = 1.1
mySim.r2 = 1.2
mySim.populationHistory = 0

mySim.in0 = 10
mySim.in1 = 0
mySim.in2 = 0



dataPointCount = 40

simsPerDataPoint = 20

#set some pulsing parameters
drug_strength=0.6
freq=1.0

#pick initial r1, to find the starting point in r1/r2 space
r1_origin=1.25
r2_origin=(1.0-mySim.u2)*r1_origin

mut2= math.pow((1.0-mySim.u2),2) #define a term that appears alot in finding the initial r1 and r2
#Get the starting r1 and r2 - using formula derived in notebook
init_r1=math.pow( ( mut2 / (1.0+mut2) ), 0.5) * drug_strength/5.0 + r1_origin
init_r2=math.pow( (1.0/(1.0+mut2) ) ,0.5)* drug_strength/2.0 + r2_origin


angle_of_closest_approach=math.atan(1-mySim.u2) #find the angle that gives the closest way to get to 
angle_of_closest_approach*=180/math.pi

#do the simulation
fixTime=[]
angleRange=numpy.linspace(0,math.pi/2.0, num=dataPointCount)
for angle in angleRange:
    #allocate to presim so that the pulse is varied
    
    mySim.preSim=multiple_pulse.multiple_pulse (mySim,freq, drug_strength, angle, init_r1, init_r2, r1_offset=0.0, r2_offset=0.0, r1_width=1.0, r2_width=1.0)
    fixTime_term=0.0    
    print("Angle in radians - supplied to the simulator={0}".format(angle))
    print("Current Data Point = {0}/{1} ({2}%)".format(angle*180/math.pi, 90, 100.0 * float(2.0*angle)/math.pi))
    for i in range(0,simsPerDataPoint):
        mySim.Simulate()
        if mySim.Fixated():
            fixTime_term += mySim.curTime
        else:
            fixTime_term += mySim.timeLimit
            print("WARNING! - system did not fixate")
    average_fix=fixTime_term/simsPerDataPoint
    fixTime.append(average_fix)

angleRange*=180/math.pi #express agnles in degrees
plt.plot(angleRange,fixTime)
plt.xlabel("Angle")
plt.ylabel("Fixation Time")

dirac_angle=[]
lower, upper = plt.ylim()
dirac_time=[lower,upper]
dirac_angle.append(angle_of_closest_approach)
dirac_angle.append(angle_of_closest_approach)

plt.plot(dirac_angle,dirac_time)
plt.show





















