# -*- coding: utf-8 -*-
"""
Created on Thu Nov 06 13:48:28 2014

@author: Connor
"""

import TwoSpecies
import math
import anfixtime

import SimTools

import matplotlib.pyplot as plt

sim3 = SimTools.Gillespie(10)
sim3.timeLimit = 10000
sim3.u1 = 0.1
sim3.u2 = 0.01
sim3.in0 = 10

sim3.r0 = 1.5
sim3.r1 = 1.0
sim3.r2 = 0.8

sim3.n1DidFix = 0
sim3.n1FixTime = 0.0

def RecordFix1(sim):
    if sim.n1DidFix == 0 and sim.n1 == sim.N:
        sim.n1FixTime = sim.curTime
        sim.n1DidFix = 1

sim3.preSim = RecordFix1

dataPointCount = 30

minu1 = 0.01
maxu1 = 0.2

sdp = 1000 #simsperdatapoint

dataX = []
dataFix = []
dataFix1 = []

dataTheory = []

#This simulator considers the transition from type0 to type1
sim01 = TwoSpecies.Gillespie(10)
sim01.u1 = sim3.u1

sim01.r0 = sim3.r0
sim01.r1 = sim3.r1

sim01.ij = 0

#Now consider transition from type1 to type2
sim12 = TwoSpecies.Gillespie(10)
sim12.u1 = sim3.u2

sim12.r0 = sim3.r1
sim12.r1 = sim3.r2

sim01.ij = 0

dataTheory01 = []
dataTheory12 = []


for i in range(0,dataPointCount):
    fixTime = 0.0
    fixTime1 = 0.0
    sim3.u1 = float(maxu1 - minu1) * float(i) / (dataPointCount - 1.0) + minu1
    print(sim3.u1)
    
    for sim in range(0, sdp):
        sim3.n1DidFix = 0
        sim3.n1FixTime = 0.0
        sim3.Simulate()
        if sim3.Fixated():
            fixTime += sim3.curTime
        else:            
            print("WARNING DIDNT REACH FIXATION")
        '''if sim3.n1DidFix:
            fixTime1 += sim3.n1FixTime
        else:
            print("WARNING n1 DIDNT FIX")'''
    
    dataX.append(sim3.u1)
    dataFix.append(fixTime / sdp)
    
    dataFix1.append(fixTime1 / sdp)
    
    sim01.u1 = sim3.u1
    
    theory01 = anfixtime.GetFixTimeJ(sim01, sim01.ij)
    theory12 = anfixtime.GetFixTimeJ(sim12, sim12.ij)

    print("FIX 0->1 = {0}".format(theory01))
    print("FIX 1->2 = {0}".format(theory12))
    
    dataTheory.append(theory01 + theory12)
    dataTheory01.append(theory01)
    dataTheory12.append(theory12)

plt.plot(dataX,dataFix, '^')
#plt.plot(dataX,dataFix1, 'o')

plt.plot(dataX,dataTheory)
#plt.plot(dataX,dataTheory01, '--')
#plt.plot(dataX,dataTheory12, ':')