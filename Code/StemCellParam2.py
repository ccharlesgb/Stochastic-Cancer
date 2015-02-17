# -*- coding: utf-8 -*-
"""
Created on Thu Feb 05 13:57:23 2015

@author: Connor
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Feb 03 12:01:50 2015

@author: Connor
"""

import random

class StemCellHist:
    def __init__(self):
       self.ClearFrames()
        
    def ClearFrames(self):
        self.tHist = []
        self.n0Hist = []
        self.m0Hist = []
        self.n1Hist = []
        self.m1Hist = []
        
    def RecordFrame(self, time, param):
        if random.random() < 1.1:
            self.tHist.append(time)
            self.n0Hist.append(param.n[0])
            self.m0Hist.append(param.n[1])
            self.n1Hist.append(param.n[2])
            self.m1Hist.append(param.n[3])
        
    def GetDictionary(self):
        runDict = dict()
        runDict["time"] = self.tHist
        runDict["n0"] = self.n0Hist
        runDict["m0"] = self.m0Hist
        runDict["n1"] = self.n1Hist
        runDict["m1"] = self.m1Hist
        
        return runDict

#Define all the parameters for our stem cell model
class StemCellParam:
    def __init__(self):
        #Define Population Counts

        self.n = [0,0,0,0]
        self.n0 = [200,0,1e5,0]

        #Define stem cell reproduction rates
        self.rn = 10.0 #Normal stem cell reproduction
        self.rm = 1.0 #Cancer stem cell reproduction
        self.an = 1.065e7 #Differentiated normal cell birth rate
        self.am = 1.065e7 #Differentiated cacner cell death rate
        self.dn0 = 0.002 #Stem cell death rate healthy
        self.dm0 = 0.002 #Stem cell death rate cancer
        self.dn1 = 0.213 #Differentiated cell death rate healthy
        self.dm1 = 0.213 #Differentiated cell death rate cancerous        
        #Define Homeostasis Paramaters
        self.cn = 0.75e-3
        self.cm = 0.38e-3
        
        self.Eventn0 = [1,0,0,0]
        self.Eventn0_ = [-1,0,0,0]
        self.Eventm0 = [0,1,0,0]
        self.Eventm0_ = [0,-1,0,0]
        
    def Reset(self):
        #DO IT LIKE THIS
        for i in range(0,len(self.n)):
            self.n[i] = self.n0[i]
        
    def GetPhiNormal(self):
        return 1.0 / (1.0 + self.cn*(self.n[0] + self.n[1]))
    
    def GetPhiCancer(self):
        return 1.0 / (1.0 + self.cm*(self.n[0] + self.n[1]))
            
    #Rate at which normal stem cells reproduce
    def GetTn0(self):
        return self.n[0] * self.rn * self.GetPhiNormal()
    #Rate at which normal stem cells die
    def GetTn0_(self):
        return self.n[0] * self.dn0
        
    #Rate at which cancer stem cells reproduce
    def GetTm0(self):
        return self.n[1] * self.rm * self.GetPhiCancer()
    #Rate at which cancer stem cells die
    def GetTm0_(self):
        return self.n[1] * self.dm0
        
    def Hook(self, gillespie):
        gillespie.AddCallback(self.GetTn0, self.Eventn0)
        gillespie.AddCallback(self.GetTn0_, self.Eventn0_)
        gillespie.AddCallback(self.GetTm0, self.Eventm0)
        gillespie.AddCallback(self.GetTm0_, self.Eventm0_)
        
    def PreSim(self, gillespie):
        a = 0
        
    def PostSim(self, gillespie):
        self.n[2] += gillespie.timeStep * self.n2_derivative()
        self.n[3] += gillespie.timeStep * self.n3_derivative()       

    def n2_derivative(self):
        result = self.an*self.n[0] - self.dn1*self.n[2]
        return result
    
    def n3_derivative(self):
        result = self.am*self.n[1] - self.dm1*self.n[3]
        return result
        