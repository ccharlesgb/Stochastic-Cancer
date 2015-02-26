# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 15:59:40 2015

@author: Jonny
"""
import math
import numpy as np
import random
import scipy.misc as scimisc

class BatchResult:
    def __init__(self):
        self.simCount = 0
        self.avgFixTime = 0.0
        self.Reset()

    def Reset(self):
        self.simCount = 0
        self.avgFixTime = 0.0

class wright_fisher:
    def __init__(self, cellTypes):
        self.popSize = 100

        self.r = []
        self.useApproxTheta = 0
        self.u = 0.01
        self.cellTypes = cellTypes
        self.N = []
        self.iN = []
        self.prob_vector = []
        self.d = 100
        
        for i in range(0,self.cellTypes):
            self.N.append(0)            
            self.iN.append(0)
            self.r.append(1.0)
            self.prob_vector.append(0.0)
        
        self.iN[0] = self.popSize
        
        self.printProgress = 0
        self.curStep = 0        
        self.isFixated = 0
        self.stepLimit = 2000
        self.history = 0
        
    def SetCompoundFitness(self,s):
        for i in range(0,self.cellTypes):
            self.r[i] = math.pow(1.0 + s, i)
        
    def CacheCombinations(self):
        self.combinations = [[]]
        for i in range(0,self.cellTypes):
            self.combinations.append([])
            for j in range(0,self.cellTypes):
                comb = scimisc.comb(self.d-i, j-i)
                self.combinations[i].append(comb)
        
    def reset(self):
        self.CacheCombinations()
        #self.iN = [0]*self.cellTypes
        #self.iN[0] = self.popSize 
        self.popSize = 0
        for i in range(0,self.cellTypes):
            self.N[i] = self.iN[i]
            self.popSize += self.iN[i]
        self.curStep = 0 
        self.isFixated = 0       
        self.nextProgressFrac = 0.0
    
    def GetXi(self, i):     
        result = float(self.N[i])/float(self.popSize)       
        return result
    
    def GetAvgFitness(self):
        bottom = 0.0       
        for l in range(0,self.cellTypes):
            bottom += self.r[l]*self.N[l]
        return bottom
        
    def GetThetaj(self,j):
        if self.useApproxTheta == 1:
            avgFit = 0.0
            for l in range(0, self.cellTypes):
                avgFit += self.r[l]*self.N[l]
            #print(avgFit)
            theta_j = (self.r[j] * self.N[j])/avgFit
            if j > 0:
                theta_j += float(self.u * (self.cellTypes - j + 1) * self.r[j-1]*self.N[j-1] )/ avgFit
                
            return theta_j
        else:
            summation = 0.0
            avgFit = self.GetAvgFitness()
            for i in range(0, j+1):
                #summation += scimisc.comb(self.cellTypes - i , j - i )*math.pow(self.u, j-i)*math.pow(1-self.u, self.cellTypes-j)*self.GetFitnessRatio(i)
                summation += self.combinations[i][j]*math.pow(self.u, j-i)*math.pow(1-self.u, self.d-j)*(self.r[i] * self.N[i])/avgFit
            return summation
        
    def UpdateProbVector(self):    
        probSum = 0.0
        for i in range(0,self.cellTypes - 1):
            self.prob_vector[i] = (self.GetThetaj(i))
            probSum += self.prob_vector[i]
        self.prob_vector[self.cellTypes - 1] = 1.0 - probSum
        
        #print(self.prob_vector)
        #normalisation = sum(self.prob_vector)
        #for i in range(0,self.cellTypes):
            #self.prob_vector[i] /= normalisation
         
        for i in range(0,self.cellTypes):        
            if(self.prob_vector[i] == 1.0):
                    self.isFixated = 1      
                    print("The system is fixed. Took {0} steps".format(self.curStep))
    
    def SetHistory(self, hist):
        self.history = hist
    
    
    def Simulate(self):  
        self.reset()
        if(self.history != 0):        
            self.history.RecordFrame(self)
        
        while(self.curStep < self.stepLimit and self.isFixated != 1):
            if (self.printProgress and float(self.curStep) / self.stepLimit > self.nextProgressFrac):
                print(float(self.curStep) / self.stepLimit * 100.0)
                self.nextProgressFrac += 0.1
            self.UpdateProbVector()
            self.N = np.random.multinomial(self.popSize, self.prob_vector)           
            self.curStep += 1
            
            if self.N[self.cellTypes-1] >= 1:
                self.isFixated = 1
            
            if(self.history != 0):        
                self.history.RecordFrame(self)
            
    def SimulateBatch(self, simCount):
        res = BatchResult()
        res.simCount = simCount
        for i in range(0, simCount):
            self.Simulate()
            res.avgFixTime += self.curStep
            
        res.avgFixTime = float(res.avgFixTime) / simCount
        return res
            
    def AnalyticalWaitingTime(self):
        return 0.0
    
class wf_hist:
    def __init__(self, cellTypes):
        self.cellTypes = cellTypes
        self.ClearFrames()
       
    def ClearFrames(self):
        self.histArray = dict()
        self.stepHist = []
        self.thetajHist = dict()
        self.avgJHist = []
        for i in range(0, self.cellTypes):
            self.histArray[i] = []
            self.thetajHist[i] = []
            
    def RecordFrame(self, sim):
        if random.random() < 1.1:
            self.stepHist.append(sim.curStep)
            totalJ = 0.0
            for i in range(0, self.cellTypes):
                self.histArray[i].append(sim.N[i])
                self.thetajHist[i].append(sim.prob_vector[i])
                totalJ += i * float(sim.N[i]/sim.popSize)
            self.avgJHist.append(totalJ)
         
    def GetDictionary(self):
        runDict = dict()
        return runDict
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        