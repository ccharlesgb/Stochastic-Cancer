# -*- coding: utf-8 -*-
"""
Created on Tue Feb 03 13:10:36 2015

@author: Connor
"""

import math

class CarcinogenNHist:
    def __init__(self, typeCount):
        self.typeCount = typeCount
        self.ClearFrames()
       
    def ClearFrames(self):
        self.tHist = []
        self.histArray = dict()
        for i in range(0, self.typeCount):
            self.histArray[i] = []

    def RecordFrame(self, time, param):
            self.tHist.append(time)
            for i in range(0, self.typeCount):
                self.histArray[i].append(param.n[i])
         
    def GetDictionary(self):
        runDict = dict()
        return runDict

#Define all the parameters for our stem cell model
class CarcinogenNParam:
    def __init__(self, typeCount):
        self.typeCount = typeCount
        self.n = []
        self.n0 = []
        self.r = []
        self.u = []
        self.d = typeCount
        
        for i in range(0, typeCount):
            self.n.append(0)
            self.n0.append(0)
            self.r.append(1.0)
            self.u.append(0.0001)
            
        self.n[0] = 10
        self.avgFit = 0.0
        
    def SetTypeCount(typeCount):
        return
        
    def Reset(self):
        self.N = 0        
        for pop in range(0,len(self.n)):        
            self.n[pop] = self.n0[pop]
            self.N += self.n0[pop]
            
    def EventTIJ(self, i,j):
        arr = []
        for x in range(0,self.typeCount):
            arr.append(0)
            if x == i:
                arr[x] = -1
            elif x == j:
                arr[x] = 1
        return arr
    
    def Hook(self, gillespie):
        for i in range(0,self.typeCount):
            for j in range(0, self.typeCount):
                if i != j:
                    gillespie.AddCallback(0, self.EventTIJ(i,j))

    #Todo calculate this once
    def GetAvgFit(self):
        tot = 0.0
        for i in range(0, self.typeCount):
            tot += self.r[i] * self.n[i]
        return tot
    
    #Reaction probability for cell from 1->0
    def GetTIJ(self, i, j):
        if self.n[i] == 0: #Quick get out case
            return 0.0
        top = 0.0
        u_j = self.u[j] * (self.d - j + 1)
        u_jm1 = self.u[j-1] * (self.d - (j-1) + 1)
        if j == 0:
            top = self.n[i] * (self.r[j] * (1 - u_j) * self.n[j])
        elif j == self.typeCount - 1:
            top = self.n[i] * (self.r[j]*self.n[j] + self.r[j-1] * u_jm1 * self.n[j-1])
        else:
            top = self.n[i] * (self.r[j]*(1 - u_j)*self.n[j] + self.r[j-1]*u_jm1*self.n[j-1])
        rate = top / self.avgFit
        #if rate > 1e6:
            #print("Large Rate T{0}>{1}: {2} Avg Fit: {3} n_{4} = {5} n_{6} = {7}".format(i,j,rate, self.avgFit, i, self.n[i],j,self.n[j]))
        return rate
        
    def PreSim(self, gillespie):
        self.avgFit = self.GetAvgFit()
    
    def PostSim(self, gillespie):
        return
