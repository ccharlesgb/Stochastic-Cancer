# -*- coding: utf-8 -*-
"""
Created on Tue Feb 03 13:10:36 2015

@author: Connor
"""

class CarcinogenHist:
    def __init__(self):
       self.ClearFrames()
        
    def ClearFrames(self):
        self.tHist = []
        self.n0Hist = []
        self.n1Hist = []
        self.n2Hist = []
        
    def RecordFrame(self, time, param):
        self.tHist.append(time)
        self.n0Hist.append(param.n[0])
        self.n1Hist.append(param.n[1])
        self.n2Hist.append(param.n[2])
        
    def GetDictionary(self):
        runDict = dict()
        runDict["time"] = self.tHist
        runDict["n0"] = self.n0Hist
        runDict["n1"] = self.n1Hist
        runDict["n2"] = self.n2Hist
        
        return runDict

#Define all the parameters for our stem cell model
class CarcinogenParam:
    def __init__(self):
        self.n = [0,0,0]
        
        self.n0 = [10,0,0]
        
        self.r0 = 1.0
        self.r1 = 1.0
        self.r2 = 1.0
        
        self.u1 = 0.0
        self.u2 = 0.0
        
        self.addRate = 0.01
        
        self.c0 = 10.0
        self.c1 = 10.0
        self.c2 = 10.0
        
        self.EventT10 = [+1,-1,0]    
        self.EventT01 = [-1,+1,0]
        self.EventT20 = [+1,0,-1]
        self.EventT02 = [-1,0,+1]
        self.EventT12 = [0,-1,+1]
        self.EventT21 = [0,+1,-1]
        
        self.EventT00 = [1,0,0]
        self.EventT11 = [0,1,0]
        self.EventT22 = [0,0,1]
        
    def AvgFitness(self):
        return (self.r0*self.n[0] + self.r1*self.n[1] + self.r2 * self.n[2])
        
    def Reset(self):
        self.N = 0        
        for pop in range(0,len(self.n)):        
            self.n[pop] = self.n0[pop]
            self.N = self.n0[0] + self.n0[1] + self.n0[2]
    
    def Hook(self, gillespie):
        gillespie.AddCallback(self.GetT10, self.EventT10)
        gillespie.AddCallback(self.GetT01, self.EventT01)
        gillespie.AddCallback(self.GetT20, self.EventT20)
        gillespie.AddCallback(self.GetT02, self.EventT02)
        gillespie.AddCallback(self.GetT12, self.EventT12)
        gillespie.AddCallback(self.GetT21, self.EventT21)
        
        gillespie.AddCallback(self.GetT00, self.EventT00)
        gillespie.AddCallback(self.GetT11, self.EventT11)
        gillespie.AddCallback(self.GetT22, self.EventT22)
    
    #Reaction probability for cell from 1->0
    def GetT10(self):
        top = (1.0 - self.u1) * self.r0 * float(self.n[0]) * self.n[1]
        return top / self.AvgFitness()
    
    def GetT20(self):
        top = (1.0 - self.u1) * self.r0 * float(self.n[0]) * self.n[2]
        return top / self.AvgFitness()
    
    def GetT01(self):
        top = (self.u1 * self.r0 * self.n[0] + (1.0-self.u2) * self.r1 * self.n[1]) * self.n[0]
        return top / self.AvgFitness()
    
    def GetT21(self):
        top = (self.u1 * self.r0 * self.n[0] + (1.0-self.u2) * self.r1 * self.n[1]) * self.n[2]
        return top / self.AvgFitness()
    
    def GetT02(self):
        top = (self.u2 * self.r1 * self.n[1] + self.r2 * self.n[2]) * self.n[0]
        return top / self.AvgFitness()
        
    def GetT12(self):
        top = (self.u2 * self.r1 * self.n[1] + self.r2 * self.n[2]) * self.n[1]
        return top / self.AvgFitness()
    
    def GetT00(self):
        top = self.n[0] * self.r0 * self.addRate
        return top / (1.0 + self.c0 * self.N * self.N)
        
    def GetT11(self):
        top = self.n[1] * self.r1 * self.addRate
        return top / (1.0 + self.c1 * self.N * self.N)
        
    def GetT22(self):
        top = self.n[2] * self.r2 * self.addRate
        return top / (1.0 + self.c2 * self.N * self.N)
        
    def PreSim(self, gillespie):
        return        
    
    def PostSim(self, gillespie):
        return
