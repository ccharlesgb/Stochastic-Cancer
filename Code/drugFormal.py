# -*- coding: utf-8 -*-
"""
Created on Mon Nov 03 16:43:17 2014

@author: Connor
"""

#Drug Testing
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import SimTools
import time
import math
import numpy

def FixedPointExists(r1, r2, u2):
    condition = r2 / (1.0 - u2)
    return r1 > condition and r2 < condition
    
def BoundaryFP(r1, r2, u2):
    exists = (r1 > r2/(1.0-u2)) and (r2 < r2/(1.0-u2))
    if exists == 0:
        fp = []
        fp.append(numpy.nan)
        fp.append(numpy.nan)
        fp.append(numpy.nan)
        return fp
    fp = []
    fp.append(numpy.nan) #Append n0
    fp.append(1.0 - u2*r1/(r1-r2))
    fp.append(u2*r1/(r1-r2))
    return fp
    
def ReactiveFP(r0, r1, r2, u1, u2):
    denom = u2*r1*(r0-r2) + (r0-r1)*((1.0-u1)*r0 - r2)
    exists = ((1.0-u1)*r0 > (1.0-u2)*r1) and ((1.0-u1)*r0 > r2)
    if exists == 0:
        fp = []
        fp.append(numpy.nan)
        fp.append(numpy.nan)
        fp.append(numpy.nan)
        return fp
    numer1 = ((1.0-u1)*r0 + u2*r1 - r2)*u1*r0
    fp = []
    fp.append(1.0 - numer1/denom)
    fp.append((((1.0-u1)*r0-r2)*u1*r0)/denom)
    fp.append(u1*u2*r0*r1/denom)
    return fp
    
def FixedPointStable(r0, r1, u2):
    return (1.0 - u2)*r1 > (1.0 - u2) * r0
def FixedPointSaddle(r0,r1, u1, u2):
    exists = (1.0 - u2)*r1 < (1.0 - u1)*r0
    if exists == 0:
        return -1.0
    return 

#Initialize the Gillespie simulator with 10 cells
mySim = SimTools.Gillespie(10)

mySim.in0 = 10
mySim.timeLimit = 100.0
mySim.u1 = 0.1
mySim.u2 = 0.1
mySim.r1 = 1.1
mySim.r2 = 1.2
mySim.populationHistory = 0

mySim.in0 = 10
mySim.in1 = 0
mySim.in2 = 0

dataPointCount = 15

simsPerDataPoint = 2000

maxArea = 1.0

drugAng = 0.0
maxAng  = 90.0

#Place r1 and r2 on the boundary of the FP
r1Origin = 1.25
r2Origin = (1-mySim.u2)*r1Origin

normGrad = -1.0/(1.0 - mySim.u2)

print("NormGrad: {0}".format(normGrad))

dist = 0.1
r1Origin = r1Origin - dist
r2Origin = r2Origin + normGrad * (-dist)

print("r1Origin = {0} r2Origin = {1}".format(r1Origin, r2Origin))

for g in range(0,5):
        
    #Initialize the array with default values
    dataPointsX = []
    dataPointsY = []
    
    mySim.u2 = (0.1-0.01) * float(g)/4 + 0.01
    
    drugEffect = (0.9-0.5)*float(0)/4 + 0.5
    
    for curPoint in range(0, dataPointCount):
        startTime = time.clock() #Algorithm benchmarking
        
        fixationCounts = 0
        
        drugAng = (float(curPoint) / (dataPointCount - 1)) * (math.pi / 2.0)
        
        mySim.r1 = r1Origin + math.cos(drugAng) * drugEffect
        mySim.r2 = r2Origin - math.sin(drugAng) * drugEffect
    
        print("Current Data Point = {0}/{1} ({2}%)".format(curPoint + 1, dataPointCount, 100.0 * float(curPoint+1.0)/dataPointCount))
        print("R1 = {0}".format(mySim.r1))
        print("R2 = {0}".format(mySim.r2))
        
        #Perform many simulations to get an accurate probability of the fixation probability
        for sim in range(0, simsPerDataPoint):
            mySim.Simulate()
                
            if mySim.n2 == mySim.N: #The simulation ended with fixation
                fixationCounts += 1
        #Once the loop is done get the fraction of fixations for this r1
        dataPointsX.append(drugAng * 180.0 / math.pi)
        dataPointsY.append(float(fixationCounts) / float(simsPerDataPoint))
    
        print("Complete (Took {:.{s}f} seconds)".format(time.clock() - startTime, s=2))
    
    err = 1.0 / math.sqrt(simsPerDataPoint)
    plt.errorbar(dataPointsX, dataPointsY,yerr = err, label = "D_dose = {0}".format(drugEffect))
    plt.xlabel("Drug Angle (deg")
    plt.ylabel("Fixation")
    plt.legend()
    
plt.show()

        
