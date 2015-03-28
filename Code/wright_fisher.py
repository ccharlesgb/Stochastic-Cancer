# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 15:59:40 2015

@author: Jonny
"""
import math
import numpy as np
import random
import scipy.misc as scimisc


class wright_fisher_params:
    def __init__(self, cellTypes):
        self.r = []
        self.u = []
        self.useApproxTheta = 0
        self.cellTypes = cellTypes
        self.N = []
        self.iN = []
        self.d = 100
        self.uNotConst = 0
        
        #growing populations
        self.a = 0
        self.b = 0
                
        
        for i in range(0,self.cellTypes):
            self.N.append(0)            
            self.iN.append(0)
            self.r.append(1.0)
            self.u.append(1e-7)

        
        self.iN[0] = 100
        self.popSize = 100
    
    def UpdatePopSize(self, t):
        popSize = round(self.iN[0]* math.exp(self.b * t))
        self.popSize = popSize
             
    def PreSim(self):
        self.UpdatePopSize()
        
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
        elif self.uNotConst == 1:
            summation = 0.0
            avgFit = self.GetAvgFitness()
            for i in range(0, j+1):
                summation += self.combinations[i][j]*math.pow(self.u[i], j-i)*math.pow(1-self.u[i], self.d-j)*(self.r[i] * self.N[i])/avgFit
            return summation
        else:
            summation = 0.0
            avgFit = self.GetAvgFitness()
            for i in range(0, j+1):
                #summation += scimisc.comb(self.cellTypes - i , j - i )*math.pow(self.u, j-i)*math.pow(1-self.u, self.cellTypes-j)*self.GetFitnessRatio(i)
                summation += self.combinations[i][j]*math.pow(self.u[0], j-i)*math.pow(1-self.u[0], self.d-j)*(self.r[i] * self.N[i])/avgFit
            return summation
            
        
    def CacheCombinations(self):
        self.combinations = [[]]
        for i in range(0,self.cellTypes):
            self.combinations.append([])
            for j in range(0,self.cellTypes):
                comb = scimisc.comb(self.d-i, j-i)
                self.combinations[i].append(comb)    
    
    def Reset(self):   
        self.popSize = 0
        for i in range(0,self.cellTypes):
            self.N[i] = self.iN[i]
            self.popSize += self.iN[i]
        self.CacheCombinations()
        
    def SetCompoundFitness(self,s):
        for i in range(0,self.cellTypes):
            self.r[i] = math.pow(1.0 + s, i)

    def GetXJ(self,j):
        return self.N[j] / self.popSize            
    
    def GetXJDot(self,j, xj):
        s = self.r[1] - 1.00
        avgJ = 0.0
        for i in range(0, self.cellTypes):
            avgJ += i * xj[i]
        deriv = 0.0
        if j == 0:
            deriv = self.u[0] * (self.d - j)*xj[j]
        else:
            deriv = self.u[0] * ((self.d - j + 1) * xj[j-1] - (self.d - j)*xj[j])
        deriv += s * xj[j]*(j - avgJ)
        return deriv

    def AnalyticalWaitingTime(self):
        s = self.r[1] - 1.00
        numer = (self.cellTypes-1.0) * math.pow(math.log(s/(self.u[0]*self.d)),2.0)
        denom = 2.0 * s * math.log(self.popSize)
        return float(numer)/denom



class BatchResult:
    def __init__(self):
        self.simCount = 0
        self.avgFixTime = 0.0
        self.Reset()

    def Reset(self):
        self.simCount = 0
        self.avgFixTime = 0.0
    

class wright_fisher:
    def __init__(self):
        self.printProgress = 0
        self.curStep = 0        
        self.isFixated = 0
        self.stepLimit = 2000
        self.history = 0
        self.params = 0
        self.prob_vector = []
        self.nextBatchProgressFrac = 0.1
        self.stopAtAppear = 1
        self.reset()
        
    def reset(self):
        if self.params != 0:
            self.params.Reset()
        
            for i in range(0,self.params.cellTypes):
                self.prob_vector.append(0.0)
                
        self.curStep = 0 
        self.isFixated = 0       
        self.nextProgressFrac = 0.1
        if self.history != 0:
            self.history.ClearFrames()
    
    def GetXi(self, i):     
        result = float(self.N[i])/float(self.popSize)       
        return result
        
    def UpdateProbVector(self):    
        probSum = 0.0
        for i in range(0,self.params.cellTypes - 1):
            self.prob_vector[i] = (self.params.GetThetaj(i))
            probSum += self.prob_vector[i]
        self.prob_vector[self.params.cellTypes - 1] = 1.0 - probSum
        
        #print(self.prob_vector)
        #normalisation = sum(self.prob_vector)
        #for i in range(0,self.cellTypes):
            #self.prob_vector[i] /= normalisation
         
        for i in range(0,self.params.cellTypes):        
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
            if (self.printProgress >= 2 and float(self.curStep) / self.stepLimit > self.nextProgressFrac):
                print(int(float(self.curStep) / self.stepLimit * 100.0)),
                self.nextProgressFrac += 0.1
            if (self.params != 0):
                self.params.PreSim()
            self.UpdateProbVector()
            self.params.N = np.random.multinomial(self.params.popSize, self.prob_vector)
            self.curStep += 1
            
            if self.stopAtAppear == 1 and self.params.N[self.params.cellTypes-1] >= 1:
                self.isFixated = 1
            
            if(self.history != 0):
                self.history.RecordFrame(self)
        
        if self.printProgress >= 2:
            print('100')

    def SimulateBatch(self, simCount):
        res = BatchResult()
        res.simCount = simCount
        if simCount <= 0:
            return res
        
        print("Batch %: 0"),
        for i in range(0, simCount):
            #Print Batch progress
            if (self.printProgress >= 1 and float(i) / simCount > self.nextBatchProgressFrac):
                print(int(float(i) / simCount * 100.0)),
                self.nextBatchProgressFrac += 0.1
                
            self.Simulate()
            res.avgFixTime += self.curStep
        
        if self.printProgress >= 1:
            print('100')      
            self.nextBatchProgressFrac = 0.1
        
        res.avgFixTime = float(res.avgFixTime) / simCount
        return res

class wf_hist:
    def __init__(self, cellTypes):
        self.cellTypes = cellTypes
        self.ClearFrames()
       
    def ClearFrames(self):
        self.histArray = dict()
        self.stepHist = []
        self.yearHist = []
        self.thetajHist = dict()
        self.avgJHist = []
        self.avgSJHist = []
        for i in range(0, self.cellTypes):
            self.histArray[i] = []
            self.thetajHist[i] = []
            
    def RecordFrame(self, sim):
        self.stepHist.append(sim.curStep)
        self.yearHist.append(float(sim.curStep) / 365.0)
        totalJ = 0.0
        totalSJ = 0.0            
        for i in range(0, self.cellTypes):
            self.histArray[i].append(sim.params.N[i])
            self.thetajHist[i].append(sim.prob_vector[i])
            totalJ += i * float(sim.params.N[i]/sim.params.popSize)
            totalSJ += (sim.params.r[i] - 1.0 )* float(sim.params.N[i]/sim.params.popSize)
        self.avgJHist.append(totalJ)
        self.avgSJHist.append(totalSJ)
         
    def GetDictionary(self):
        runDict = dict()
        return runDict
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        