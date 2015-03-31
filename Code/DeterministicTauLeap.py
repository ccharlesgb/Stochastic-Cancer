# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 16:25:50 2015

@author: Connor
"""

import math

class Sim: 
    def __init__(self, typeCount):
        self.history = 0
        self.params = 0
        self.typeCount = typeCount
        self.stopAtAppear = 0
        self.n = [] #Our state vector
        self.x = [] #Concentrations
        self.x_dot = []
        for i in range(0,self.typeCount):
            self.n.append(0)
            self.x.append(0.0)
            self.x_dot.append(0.0)
        
        self.timeLimit = 1e4
        self.timeStep = 1e-1
        self.Reset()
   
        
    def Reset(self):
        self.avgJ = 0.0
        self.simSteps = 0
        self.curTime = 0.0
        
        if self.history != 0:
            self.history.ClearFrames()
            
        if self.params  != 0:
            for i in range(0,self.typeCount):
                self.n[i] = self.params.n0[i]
                self.x[i] = float(self.n[i]) / self.params.N
    
    def GetXDot(self, i):
        u_minus = 0.0
        x_minus = 0.0
        s = self.params.r[1] - self.params.r[0]
        d = self.params.d
        if i > 0:
            x_minus = self.x[i-1]
            u_minus = self.params.u[i-1]

        u = self.params.u[i]
        if (i+1) >= self.typeCount:
            u = 0.0
        mut_term = u_minus * ((d-i+1)*x_minus) - u * (d-i)*self.x[i]
        
        repo_term = s * self.x[i] * (i - self.avgJ)
        
        return mut_term + repo_term
    
    def UpdateAvgJ(self):
        self.avgJ = 0.0
        for i in range(0, self.typeCount):
            self.avgJ = self.avgJ + (i * self.x[i])
            
    def UpdateStateVector(self):
        for i in range(0, self.typeCount):
            self.n[i] = math.ceil(self.x[i] * self.params.N)
    
    def Integrate(self):
        self.Reset()
        steps = int(self.timeLimit / self.timeStep)
        
        for st in range(0, steps):
            self.UpdateAvgJ()
            #Update derivatives
            for i in range(0, self.typeCount):
                self.x_dot[i] = self.GetXDot(i)
            for i in range(0, self.typeCount):
                self.x[i] = self.x[i] + self.GetXDot(i) * self.timeStep
        
            self.UpdateStateVector()
            
            self.curTime = float(st)/steps * self.timeLimit
            if self.history != 0:
                self.history.RecordFrame(self, self.params)
            
            if self.stopAtAppear == 1:
                if self.n[self.typeCount - 1] >= 1.0:
                    break;
            else:
                if self.n[self.typeCount-1] == self.params.N: #Fixated
                   break;
                
                
        