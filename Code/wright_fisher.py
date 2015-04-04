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
    def __init__(self, typeCount):
        self.typeCount = typeCount
        self.printProgress = 0
        self.curTime = 0
        self.isFixated = 0
        self.timeLimit= 1e4
        self.history = 0
        self.params = 0
        self.prob_vector = []
        self.nextBatchProgressFrac = 0.1
        self.stopAtAppear = 1
        self.n = [] #Our state vector
        for i in range(0,self.typeCount):
            self.n.append(0)
            self.prob_vector.append(0.0)
            
        self.reset()
        
    def reset(self):
        if self.params != 0:
            self.params.Reset()
            for i in range(0,self.typeCount):
                self.prob_vector[i] = 0.0
            
            for i in range(0,self.typeCount):
                self.n[i] = self.params.n0[i]        
        
        self.curTime = 0
        self.isFixated = 0       
        self.nextProgressFrac = 0.1
        if self.history != 0:
            self.history.ClearFrames()
        
    def UpdateProbVector(self):    
        probSum = 0.0
        for i in range(0,self.typeCount - 1):
            self.prob_vector[i] = (self.params.GetThetaj(i, self.n))
            probSum += self.prob_vector[i]
        self.prob_vector[self.typeCount - 1] = 1.0 - probSum
         
        for i in range(0,self.typeCount):        
            if(self.prob_vector[i] == 1.0):
                    self.isFixated = 1      
                    print("The system is fixed. Took {0} steps".format(self.curTime))
    
    def SetHistory(self, hist):
        self.history = hist
    
    def Simulate(self):
        self.reset()
        if(self.history != 0):
            self.history.RecordFrame(self)
        #print("INITIL CURTIME", self.curTime)
        #print("INITIAL n", self.n)
        #print("INITIAL popsize",self.params.N)
        while(self.curTime < self.timeLimit and self.isFixated != 1):
            if (self.printProgress >= 2 and float(self.curTime) / self.timeLimit > self.nextProgressFrac):
                print(int(float(self.curTime) / self.timeLimit * 100.0)),
                self.nextProgressFrac += 0.1
            
            if self.params != 0:
                self.params.PreSim(self)
            
            self.UpdateProbVector()
            self.n = np.random.multinomial(self.params.N, self.prob_vector)
            self.curTime += 1
            
            if self.stopAtAppear == 1 and self.n[self.params.typeCount-1] >= 1:
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
        
        if self.printProgress >= 1:
            print("Batch %: 0"),
        for i in range(0, simCount):
            #Print Batch progress
            if (self.printProgress >= 1 and float(i) / simCount > self.nextBatchProgressFrac):
                print(int(float(i) / simCount * 100.0)),
                self.nextBatchProgressFrac += 0.1
                
            self.Simulate()
            print(self.curTime)
            res.avgFixTime += self.curTime
        
        if self.printProgress >= 1:
            print('100')      
            self.nextBatchProgressFrac = 0.1
        
        res.avgFixTime = float(res.avgFixTime) / simCount
        return res


         
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        