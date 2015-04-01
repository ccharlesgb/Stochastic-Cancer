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
        self.histArray = dict()
        self.tHist = []
        self.yearHist = []
        self.thetajHist = dict()
        self.avgJHist = []
        self.avgSJHist = []
        for i in range(0, self.typeCount):
            self.histArray[i] = []
            self.thetajHist[i] = []

    def RecordFrame(self, sim):
        self.tHist.append(sim.curTime)
        self.yearHist.append(float(sim.curTime) / 365.0)
        totalJ = 0.0
        totalSJ = 0.0            
        for i in range(0, self.typeCount):
            self.histArray[i].append(sim.n[i])
            if sim.prob_vector != 0:
                self.thetajHist[i].append(sim.prob_vector[i])
            totalJ += i * float(sim.n[i]/sim.params.N)
            totalSJ += (sim.params.r[i] - 1.0 )* float(sim.n[i]/sim.params.N)
        self.avgJHist.append(totalJ)
        self.avgSJHist.append(totalSJ)
         
    def GetDictionary(self):
        runDict = dict()
        return runDict

#Define all the parameters for our stem cell model
class Params:
    def __init__(self, typeCount):
        self.SetTypeCount(typeCount)
        self.d = 100
        self.USE_D = True
        self.avgFit = 0.0
        self.uNotConst = 0
        self.Reset()
        
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
        self.n0[0] = 1e2
        
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
    
    #Reaction probability for cell from 1->0
    def GetTIJ(self, i, j, n):
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
                top = (self.r[j] * (1.0 - u_j) * n[j])
            elif j == self.typeCount - 1:
                top = (self.r[j]*n[j] + self.r[j-1] * u_jm1 * n[j-1])
            else:
                top = (self.r[j]*(1.0 - u_j)*n[j] + self.r[j-1]*u_jm1*n[j-1])
            rate = top / self.avgFit
            self.thetaJCache[j] = rate
    
    def SetCompoundFitness(self,s):
        for i in range(0,self.typeCount):
            self.r[i] = math.pow(1.0 + s, i)    
    
    def CacheCombinations(self):
        self.combinations = [[]]
        for i in range(0,self.typeCount):
            self.combinations.append([])
            for j in range(0,self.typeCount):
                comb = scimisc.comb(self.d-i, j-i)
                self.combinations[i].append(comb)        
    
    def GetThetaj(self,j, n):
        if self.uNotConst == 1:
            summation = 0.0
            avgFit = self.avgFit
            for i in range(0, j+1):
                summation += self.combinations[i][j]*math.pow(self.u[i], j-i)*math.pow(1.0-self.u[i], self.d-j)*(self.r[i] * n[i])/avgFit
            return summation
        else:
            summation = 0.0
            avgFit = self.avgFit
            for i in range(0, j+1):
                summation += self.combinations[i][j]*math.pow(self.u[0], j-i)*math.pow(1.0-self.u[0], self.d-j)*(self.r[i] * n[i])/avgFit
            return summation
    
    def PreSim(self, gillespie):
        self.avgFit = self.GetAvgFit(gillespie.n)
        self.CacheThetaJ(gillespie.n)
    
    def PostSim(self, gillespie):
        return
