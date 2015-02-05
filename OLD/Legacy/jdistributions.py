# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 14:11:47 2014

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
N=100
mySim = TwoSpecies.Gillespie(N)
mySim.timeLimit = 3000
dataPointCount = 50
mySim.r0=1.0
mySim.r1=1.0
mySim.u1=0.06

dataPointCount = 20
simsPerDataPoint = 10000

max_time=1500
min_time=150



indiv_distrib=[]
total_distrib=[]
Time_Limit=[]
meanj=[]
for curPoint in range(0, dataPointCount):
    startTime = time.clock() #Algorithm benchmarking
    j_count=0
    mySim.timeLimit=curPoint/(dataPointCount - 1.0)*(max_time - min_time) + min_time    
    
    #successful_fixations =  0 #to exclude simulations that fix at j=0 
    #mySim.ij = int(float(curPoint)/(dataPointCount-1) * mySim.N)
    mySim.ij=0
    print("Current Data Point = {0}/{1} ({2}%)".format(curPoint + 1, dataPointCount, 100.0 * float(curPoint+1.0)/dataPointCount))
    
    
    #Perform many simulations to get an accurate probability of the fixation probabilities
    #while successful_fixations < simsPerDataPoint:    
    for i in range(0, simsPerDataPoint):
        mySim.Simulate()
            
        if mySim.j == mySim.N: #or mySim.j==0: #The simulation ended with fixation
            j_count+=mySim.j
            
        indiv_distrib.append(mySim.j)
        
    total_distrib.append([])
    total_distrib[curPoint] = copy.copy(indiv_distrib) #add each data points fix time distribution to an array containing all distribs
        
    
    #Once the loop is done get the fraction of fixations for this r1
    Time_Limit.append(curPoint/(dataPointCount - 1.0)*(max_time - min_time) + min_time)
    meanj.append(float(j_count) / (float(simsPerDataPoint))) #fixation time)
    print("Complete (Took {:.{s}f} seconds)".format(time.clock() - startTime, s=2))
    
    
for i in range(0, dataPointCount):
    hist, bins = np.histogram(total_distrib[i], bins=50) #create histogram for each data points distribution with set number of bins
    width = (bins[1] - bins[0]) 
    center = (bins[:-1] + bins[1:]) / 2
    
    plt.bar(center, hist, align='center', width=width)
    plt.title("Timelimit:{0} Average Value:{1}".format(i/(dataPointCount - 1.0)*(max_time - min_time) + min_time,sum(total_distrib[i])/len(total_distrib[i])))    
    plt.xlabel("J value")
    plt.ylabel("Counts")
    plt.show()