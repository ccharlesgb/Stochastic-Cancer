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
        
    def RecordFrame(self, time, param):
        if random.random() < 1.1:
            self.tHist.append(time)
            self.n0Hist.append(param.n0)
            self.m0Hist.append(param.m0)
        
    def GetDictionary(self):
        runDict = dict()
        runDict["time"] = self.tHist
        runDict["n0"] = self.n0Hist
        runDict["m0"] = self.m0Hist
        
        return runDict

#Define all the parameters for our stem cell model
class StemCellParam:
    def __init__(self):
        #Define Population Counts
        self.n0 = 0 #Number of HEALTHY stem cells
        self.m0 = 0 #Number of CANCER stem cells
        #Define initial conditons
        self.in0 = 200
        self.im0 = 1
        #Define stem cell reproduction rates
        self.rn = 10.0 #Normal stem cell reproduction
        self.rm = 1.0 #Cancer stem cell reproduction
        self.dn0 = 0.002 #Stem cell death rate healthy
        self.dm0 = 0.002 #Stem cell death rate cancer
        #Define Homeostasis Paramaters
        self.cn = 0.75e-3
        self.cm = 0.38e-3
        
    def Reset(self):
        self.n0 = self.in0
        self.m0 = self.im0

    def GetPhiNormal(self):
        return 1.0 / (1.0 + self.cn*(self.n0 + self.m0))
    
    def GetPhiCancer(self):
        return 1.0 / (1.0 + self.cm*(self.n0 + self.m0))
            
    #Rate at which normal stem cells reproduce
    def GetTn0(self):
        return self.n0 * self.rn * self.GetPhiNormal()
    #Rate at which normal stem cells die
    def GetTn0_(self):
        return self.n0 * self.dn0
        
    #Rate at which cancer stem cells reproduce
    def GetTm0(self):
        return self.m0 * self.rm * self.GetPhiCancer()
    #Rate at which cancer stem cells die
    def GetTm0_(self):
        return self.m0 * self.dm0
        
    #Health Cell Events
    def Eventn0(self):
        self.n0 += 1
    def Eventn0_(self):
        self.n0 = self.n0 -1
    
    #Cancer Cell Events
    def Eventm0(self):
        self.m0 += 1
    def Eventm0_(self):
        self.m0 = self.m0 - 1
        
    def Hook(self, gillespie):
        gillespie.AddCallback(self.GetTn0, self.Eventn0)
        gillespie.AddCallback(self.GetTn0_, self.Eventn0_)
        gillespie.AddCallback(self.GetTm0, self.Eventm0)
        gillespie.AddCallback(self.GetTm0_, self.Eventm0_)
    
        