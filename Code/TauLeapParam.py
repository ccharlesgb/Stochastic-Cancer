# -*- coding: utf-8 -*-
"""
Created on Mon Mar 09 14:41:29 2015

@author: Connor
"""
import math
import scipy.misc as scimisc

class Hist:
    def __init__(self, typeCount):
        self.typeCount = typeCount
        self.ClearFrames()
       
    def ClearFrames(self):
        self.tHist = []
        self.histArray = dict()
        for i in range(0, self.typeCount):
            self.histArray[i] = []

    def RecordFrame(self, sim, param):
            self.tHist.append(sim.curTime)
            for i in range(0, param.typeCount):
                self.histArray[i].append(sim.n[i])
         
    def GetDictionary(self):
        runDict = dict()
        return runDict

#Define all the parameters for our stem cell model
class Params:
    def __init__(self, typeCount):
        self.SetTypeCount(typeCount)
        self.d = 2
        self.USE_D = True
        self.avgFit = 0.0
        
    def SetTypeCount(self, typeCount):
        self.typeCount = typeCount
        self.n0 = []
        self.r = []
        self.u = []
        self.thetaJCache = []
        
        for i in range(0, typeCount):
            self.n0.append(0)
            self.r.append(1.0)
            self.u.append(0.1)
            self.thetaJCache.append(0.0)
        
    def Reset(self):
        self.N = 0        
        for pop in range(self.typeCount):
            self.N += self.n0[pop]
        self.CacheCombinations()
            
    def EventTIJ(self, i,j):
        arr = []
        for x in range(0,self.typeCount):
            arr.append(0)
            if x == i:
                arr[x] = -1
            elif x == j:
                arr[x] = 1
        return arr

    #Todo calculate this once
    def GetAvgFit(self, n):
        tot = 0.0
        for i in range(0, self.typeCount):
            tot += self.r[i] * n[i]
        return tot
    
    def Hook(self, sim):
        for i in range(0,self.typeCount):
            for j in range(0, self.typeCount):
                if i != j:
                    sim.AddStateChange(i, j, self.EventTIJ(i,j))
                    #print("Adding rate {0},{1}".format(i,j))
    
    def CacheCombinations(self):
        self.combinations = [[]]
        for i in range(0,self.typeCount):
            self.combinations.append([])
            for j in range(0,self.typeCount):
                comb = scimisc.comb(self.d-i, j-i)
                self.combinations[i].append(comb)      
    
    #Reaction probability for cell from 1->0
    def GetTIJ(self, i, j, n):
        #if rate > 1e6:
            #print("Large Rate T{0}>{1}: {2} Avg Fit: {3} n_{4} = {5} n_{6} = {7}".format(i,j,rate, self.avgFit, i, self.n[i],j,self.n[j]))
        return self.thetaJCache[j] * n[i]

    def CacheThetaJ(self, n):
        for j in range(0, self.typeCount):
            if n[j] == 0 and (j > 0 and n[j-1] == 0):
                self.thetaJCache[j] = 0.0
                continue
            top = 0.0
            u_j = self.u[j] * (self.d - j)
            u_jm1 = self.u[j-1] * (self.d - (j-1))
            if self.USE_D == False:
                u_j = self.u[j]
                u_jm1 = self.u[j-1]
            if j == 0:
                top = (self.r[j] * (1 - u_j) * n[j])
            elif j == self.typeCount - 1:
                top = (self.r[j]*n[j] + self.r[j-1] * u_jm1 * n[j-1])
            else:
                top = (self.r[j]*(1 - u_j)*n[j] + self.r[j-1]*u_jm1*n[j-1])
            rate = top / self.avgFit
            self.thetaJCache[j] = rate
    
    def PreSim(self, gillespie):
        self.avgFit = self.GetAvgFit(gillespie.n)
        self.CacheThetaJ(gillespie.n)
    
    def PostSim(self, gillespie):
        return
