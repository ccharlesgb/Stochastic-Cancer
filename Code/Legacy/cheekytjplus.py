# -*- coding: utf-8 -*-
"""
Created on Tue Oct 28 16:43:28 2014

@author: Jonny
"""
import TwoSpecies

mySim = TwoSpecies.Gillespie(10)
mySim.timeLimit = 1000000
dataPointCount = 25
#mySim.j = int(mySim.N/2)
mySim.r0=1.1
mySim.r1=1.1

def SumTJinv(sim):
    sum1=0    
    for i in range(1,sim.N):
        sim.j=i        
        sum1+=(sim.N-i)/sim.GetTJplus
    return sum1