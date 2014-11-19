# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 12:38:52 2014

@author: Connor
"""

import TwoSpecies
import math
import matplotlib.pyplot as plt

sim = TwoSpecies.Gillespie(10)
sim.ij = 0
sim.N = 1000
sim.r0 = 1.1
sim.r1 = 1.0
sim.u1 = 0.1
sim.timeLimit = 1000000.0

def FPExist(sim):
    cond = sim.r1 / (1.0 - sim.u1)
    cond2 = sim.r0 > sim.r1
    return (sim.r0 >= cond) and cond2
    
def FPJ(sim):
    j_star = -1.0
    if FPExist(sim):   
        j_star = sim.N * sim.u1 * sim.r0  / (sim.r0 - sim.r1)
    return j_star
    
dataX = []
dataY = []

dataFP = []

minu1 = 0.01
maxu1 = 0.1

dataCount = 10
sdp = 100

sim.populationHistory = 1

sim.u1 = 0.05

sim.Simulate()
plt.plot(sim.tHist, sim.jHist)

FPX = []
FPX.append(0.0)
FPX.append(sim.curTime)
FPV = []
FPV.append(FPJ(sim))
FPV.append(FPJ(sim))

plt.plot(FPX,FPV)

'''
for i in range(0, dataCount):
    sim.u1 = (maxu1 - minu1) * i / (dataCount - 1) + minu1
    
    print(sim.u1)    
    print("FIX POINT EXISTS = {0}".format(FPJ(sim)))
    fixTime = 0.0
    for isim in range(0, sdp):
        sim.Simulate()
        if sim.Fixated():
            fixTime += sim.curTime
        else:
            print("WARNING: SIM DID NOT FIXATE")

    fixTime = fixTime / sdp
    dataX.append(sim.u1)
    dataY.append(fixTime)
    
plt.plot(dataX, dataY)
     '''

