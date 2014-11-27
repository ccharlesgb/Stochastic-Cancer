# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 16:10:45 2014

@author: Jonny
"""

import systematic_transitions
import SimTools
import numpy as np
import math
import matplotlib.pyplot as plt

#y=np.zeros((10,10))
#print(y[0][0])
#initialise simulation

#general simulation parameters
mySim = SimTools.Gillespie(10)

mySim.in0 = 10
mySim.timeLimit = 4000.0
mySim.u1 = 0.1
mySim.u2 = 0.1
mySim.r0 = 1.0
mySim.r1 = 1.0
mySim.r2 = 1.0

numer = systematic_transitions.systematic(mySim)
#numer.tmax = 4000.0

DPC = 10
SPD = 1000

fixTime_sim = []

min_r1 = 0.2
max_r1 = 3.0
r1 = []
errors = []


for curPoint in range(0,DPC): #simulate some shit
    r1.append((curPoint/(DPC - 1.0))*(max_r1 - min_r1) + min_r1 )
    mySim.r1 = r1[curPoint]    
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


for curPoint in range(0,DPC):
    mySim.r1 = r1[curPoint]
    fixTime_num.append(numer.Get_fix_time(mySim) )
    
plt.errorbar(r1,fixTime_sim, errors , label = "Simulation")
plt.plot(r1,fixTime_num, label = "Numerical", marker = "o")
plt.xlabel("r1")
plt.ylabel("Fixation Time")
plt.legend()
plt.show()
    
            
            
            
            
            
            
            
            
            
            
            
            
            
    
