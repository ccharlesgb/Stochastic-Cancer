# -*- coding: utf-8 -*-
"""
Created on Sat Nov 22 15:03:45 2014

@author: Connor
"""

import SimTools
import matplotlib.pyplot as plt

sim = SimTools.Gillespie(100)
sim.r0 = 1.0
sim.r1 = 0.9
sim.r2 = 1.1

sim.in0 = 100

sim.u1 = 0.01
sim.u2 = 0.01

sim.timeLimit = 100.0

def GetAvgFit(sim):
    return xi[0] * sim.r0 + xi[1] * sim.r1 + xi[2] * sim.r2
    
def GetX0Dot(sim):
    dot = xi[0] * ((1.0 - sim.u1) * sim.r0 - GetAvgFit(sim))
    
    return dot / GetAvgFit(sim)
    
def GetX1Dot(sim):
    dot = sim.u1*sim.r0*xi[0] + ((1.0 - sim.u2)*sim.r1 - GetAvgFit(sim))*xi[1]
    
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

def Integrate(tHist, x0Hist, x1Hist, x2Hist, sim):
    sim.ResetSim()
    curT = 0.0
    x2Fixed = 0    
    
    for it in range(0, steps):
    
        if curT < smokeTime:
            sim.u1 = smoke_u1
            sim.u2 = smoke_u1
        else:
            sim.u1 = quit_u1
            sim.u2 = quit_u1
        
        #Get rate equations
        x0dot = GetX0Dot(sim)
        x1dot = GetX1Dot(sim)
        x2dot = GetX2Dot(sim)
        
        #Euler forward them
        newx0  = xi[0] + x0dot * deltaT
        newx1  = xi[1]+ x1dot * deltaT
        newx2  = xi[2] + x2dot * deltaT
    
        xi[0] = newx0
        xi[1] = newx1
        xi[2] = newx2
        
        curT = (float(it + 1) * deltaT) 
        
        tHist.append(curT)
        x0Hist.append(xi[0])
        x1Hist.append(xi[1])
        x2Hist.append(xi[2])
        
        if abs(xi[2] - 1.0) < 0.01 and x2Fixed == 0:
            print("REACHED FIXATION AT: {0}".format(curT))
            x2Fixed = 1

Integrate(tHist, x0Hist, x1Hist, x2Hist, sim)

plt.figure()
plt.subplot(2,1,1)

plt.plot(tHist, x0Hist, label = 'x0', linewidth = 2.0)
plt.plot(tHist, x1Hist, ':',label = 'x1', linewidth = 2.0)
plt.plot(tHist, x2Hist, '--', label = 'x2', linewidth = 2.0)
PlotSmokeLine(smokeTime)

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