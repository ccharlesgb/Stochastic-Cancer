# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 15:59:40 2015

@author: Jonny
"""
import math
import numpy as np
import random

class wright_fisher:
    def __init__(self):
        self.popSize = 10
        self.s = 0.01
        self.u = 0.01
        self.cellTypes = 3   
        
        self.N = []        
        self.iN = [0]*self.cellTypes
        self.iN[0] = self.cellTypes
        self.prob_vector = []
        
        self.curStep = 0        
        self.isFixated = 0
        self.stepLimit = 3000
        self.history = 0
        
    def reset(self):
        self.N = self.iN
        self.curStep = 0 
        
    def GetXi(self, i):
        result = self.N[i]/self.popSize
        return result
    
    def GetFitnessRatio(self, i):
        top = math.pow(1.0+self.s, i) * self.GetXi(i)
        bottom = 0        
        for l in range(0,self.cellTypes+1):
            bottom += math.pow(1.0 +self.s, l)*self.GetXi(l)
        return top/bottom
        
    def GetThetaj(self,j):
        summation = 0        
        for i in range(0, j):        
            summation += np.choose(self.cellTypes - i , j - i )*math.pow(self.u, j-i)*math.pow(1-self.u, self.cellTypes-j)*self.GetFitnessRatio(i)
        return summation
        
    def UpdateProbVector(self):
        for i in range(0,self.cellTypes):
            self.prob_vector[i]=self.GetThetaj(i)
        self.prob_vector /= sum(self.prob_vector)
        for i in range(0,self.cellTypes):        
            if(self.prob_vector[i] == 1.0):
                    self.isFixated = 1      
                    print("The system is fixed")
    
    def SetHistory(self, hist):
        self.history = hist
    
    
    def Simulate(self):
        self.reset()
        if(self.history != 0):        
            self.history.RecordFrame(self.N, self.curStep)

        while(self.curStep < self.stepLimit and self.isFixated != 1):
            self.UpdateProbVector()
            self.N = np.random.multinomial(self.popSize, self.prob_vector)
            self.curStep += 1
            
            if(self.history != 0):        
                self.history.RecordFrame(self.N, self.curStep)
        
class wf_hist:
    def __init__(self, cellTypes):
        self.cellTypes = cellTypes
        self.ClearFrames()
       
    def ClearFrames(self):
        self.histArray = dict()
        self.stepHist = []        
        for i in range(0, self.cellTypes):
            self.histArray[i] = []
            
    def RecordFrame(self, N, curStep):
        if random.random() < 1.1:
            self.stepHist.append(curStep)
            for i in range(0, self.cellTypes):
                self.histArray[i].append(N[i])
         
    def GetDictionary(self):
        runDict = dict()

        
        return runDict
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        