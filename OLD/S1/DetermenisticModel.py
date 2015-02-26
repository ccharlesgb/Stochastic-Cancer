# -*- coding: utf-8 -*-
"""
Created on Wed Dec 17 13:39:39 2014

@author: Connor
"""

import math

class Determenistic: 
    def __init__(self):
        self.Reset()
    
    def Reset(self):
        #Initialize concentrations
        self.xi = []
        self.xi.append(1.0)
        self.xi.append(0.0)
        self.xi.append(0.0)
        self.curT = 0.0
    
    def IsFixed(self,sim):
        exists_bound = (sim.r1 > sim.r2/(1.0-sim.u2)) and (sim.r2 < sim.r2/(1.0-sim.u2))
        exists_react = ((1.0-sim.u1)*sim.r0 > (1.0-sim.u2)*sim.r1) and ((1.0-sim.u1)*sim.r0 > sim.r2)
        return exists_bound or exists_react
            
    def GetAvgFit(self,sim):
        return (self.xi[0] * sim.r0 + self.xi[1] * sim.r1 + self.xi[2] * sim.r2) / sim.N
        
    def GetX0Dot(self,sim):
        dot = self.xi[0] * ((1.0 - sim.u1) * sim.r0 - self.GetAvgFit(sim))
        
        return dot / self.GetAvgFit(sim)
        
    def GetX1Dot(self,sim):
        dot = sim.u1*sim.r0*self.xi[0] + ((1.0 - sim.u2)*sim.r1 - self.GetAvgFit(sim))*self.xi[1]
        
        return dot / self.GetAvgFit(sim)

        
    def GetX2Dot(self,sim):
        dot = sim.u2*sim.r1*self.xi[1] + (sim.r2 - self.GetAvgFit(sim))*self.xi[2]
        
        return dot / self.GetAvgFit(sim)
        
    def Integrate(self,sim):
        deltaT = 0.2
        steps = int(sim.timeLimit / deltaT)
        sim.ResetSim()
        print("SIM N0", sim.n0)
        print("SIM N", sim.N)
        self.curT = 0.0
        x2Fixed = 0    
        
        if self.IsFixed(sim) and sim.populationHistory == 0: #We are never going to reach fixation
            self.curT = sim.timeLimit
            print("IS FIXED AT r1 = {0} r2 = {1}".format(sim.r1, sim.r2))
            return
        
        for it in range(0, steps):
            if sim.preSim != 0:
                sim.preSim(sim)
            #Get rate equations
            x0dot = self.GetX0Dot(sim)
            x1dot = self.GetX1Dot(sim)
            x2dot = self.GetX2Dot(sim)
            
            #Euler forward them
            newx0  = self.xi[0] + x0dot * deltaT
            newx1  = self.xi[1] + x1dot * deltaT
            newx2  = self.xi[2] + x2dot * deltaT
        
            self.xi[0] = newx0 
            self.xi[1] = newx1
            self.xi[2] = newx2
            
            self.curT = (float(it + 1) * deltaT)
            
            if sim.populationHistory >= 1:
                sim.tHist.append(self.curT)
                sim.n0Hist.append(round(self.xi[0] * sim.N))
                sim.n1Hist.append(round(self.xi[1] * sim.N))
                sim.n2Hist.append(round(self.xi[2] * sim.N))
            
            if math.ceil(self.xi[2] * sim.N) >= sim.N and x2Fixed == 0:
                #print("REACHED FIXATION AT: {0}".format(self.curT))
                break