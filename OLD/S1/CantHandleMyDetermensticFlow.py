# -*- coding: utf-8 -*-
"""
Created on Sat Nov 22 15:03:45 2014

@author: Connor
"""

import SimTools
import matplotlib.pyplot as plt
import MatTools
import DetermenisticModel
import math

sim = SimTools.Gillespie(1e7)
sim.populationHistory = 1
sim.r0 = 1.0
sim.r1 = math.pow(1.0 + 0.01,1)
sim.r2 = math.pow(1.0 + 0.01,2)

sim.in0 = 1e7
sim.N = 1e7

sim.u1 = 1e-7
sim.u2 = 1e-7

sim.timeLimit = 100.0

deterSim = DetermenisticModel.Determenistic()

deterSim.Integrate(sim)

'''
def GetAvgFit(sim):
    return xi[0] * sim.r0 + xi[1] * sim.r1 + xi[2] * sim.r2
    
def GetX0Dot(sim):
    dot = xi[0] * ((1.0 - sim.u1) * sim.r0 - GetAvgFit(sim))
    
    return dot / GetAvgFit(sim)
    
def GetX1Dot(sim):
    dot = sim.u1*sim.r0*xi[0] + ((1.0 - sim.u2)*sim.r1 - GetAvgFit(sim))*xi[1]
    
    return dot / GetAvgFit(sim)
    
def GetX1Dot2(sim):
    dot = xi[0] * (sim.u1 * sim.r0 * xi[0] + (1.0 - sim.u2) * sim.r1 * xi[1])
    dot += xi[2] * (sim.u1 * sim.r0 * xi[0] + (1.0 - sim.u2) * sim.r1 * xi[1])
    dot -= xi[1] * ((1.0 - sim.u1) * sim.r0 * xi[0])
    dot -= xi[1] * (sim.u2 * sim.r1 * xi[1] + sim.r2 * xi[2])
    
    return dot / GetAvgFit(sim)
    
def GetX2Dot(sim):
    dot = sim.u2*sim.r1*xi[1] + (sim.r2 - GetAvgFit(sim))*xi[2]
    
    return dot / GetAvgFit(sim)

def PlotSmokeLine(smokeT):
    xVals = []
    yVals = []
    xVals.append(smokeT)
    xVals.append(smokeT)
    yVals.append(0.0)
    yVals.append(1.0)
    plt.plot(xVals, yVals, '--')

#Euler Forward these rate equations
    
deltaT = 0.1
steps = int(sim.timeLimit / deltaT)

#Initialize concentrations
xi = []
xi.append(1.0)
xi.append(0.0)
xi.append(0.0)

#Set up hist arrays
tHist = []
x0Hist = []
x1Hist = []
x2Hist = []
tHist.append(0.0)
x0Hist.append(xi[0])
x1Hist.append(xi[1])
x2Hist.append(xi[2])

#Mutation parameters
smokeArea_u1 = 0.5
quit_u1 = 0.01
smoke_u1 = quit_u1 * 10.0

mut_factor = 20.0
smoke_u1 = quit_u1 * mut_factor
smokeTime = smokeArea_u1 / (smoke_u1 - quit_u1)

print(smoke_u1)
print(smokeTime)
'''

plt.figure()
plt.subplot(2,1,1)

plt.plot(sim.tHist, sim.n0Hist, label = 'x0', linewidth = 2.0)
plt.plot(sim.tHist, sim.n1Hist, ':',label = 'x1', linewidth = 2.0)
plt.plot(sim.tHist, sim.n2Hist, '--', label = 'x2', linewidth = 2.0)
#PlotSmokeLine(smokeTime)

plt.xlabel("time t")
plt.ylabel("Concentrations")
plt.legend()

xi = []

xi.append(1.0)
xi.append(0.0)
xi.append(0.0)

#Set up hist arrays
tHist = []
x0Hist = []
x1Hist = []
x2Hist = []

tHist.append(0.0)
x0Hist.append(xi[0])
x1Hist.append(xi[1])
x2Hist.append(xi[2])
'''
mut_factor = 2.0

smoke_u1 = quit_u1 * mut_factor
smokeTime = smokeArea_u1 / (smoke_u1 - quit_u1)

sim.ResetSim()
Integrate(tHist, x0Hist, x1Hist, x2Hist, sim)

plt.subplot(2,1,2)

plt.plot(tHist, x0Hist, label = 'x0', linewidth = 2.0)
plt.plot(tHist, x1Hist, ':',label = 'x1', linewidth = 2.0)
plt.plot(tHist, x2Hist, '--', label = 'x2', linewidth = 2.0)
PlotSmokeLine(smokeTime)

plt.xlabel("time t")
plt.ylabel("Concentrations")

plt.legend()

data = dict()
data["Deter_t"] = tHist
data["Deter_x0"] = x0Hist
data["Deter_x1"] = x1Hist
data["Deter_x2"] = x2Hist

file_name = "FP_Oscillations_determ_N={0}_r1_{1}_r2_{2}".format(sim.N, sim.r1, sim.r2)
MatTools.SaveDict(file_name, data)'''