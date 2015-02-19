# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 15:59:40 2015

@author: Jonny
"""
import math

class wright_fisher:
    def __init__(self):
        self.popSize = 10
        self.s = 0.01
        self.u = 0.01
        self.cellTypes = 3   
        
        self.N = []        
        self.iN = [0]*self.popSize
        self.iN[0] = self.popSize
        
        self.curTime = 0
    
    def reset(self):
        self.N = self.iN
        self.curTime = 0        
        
    def GetXi(self, i):
        result = self.N[i]/self.popSize
        return result
    
    def GetFitnessRatio(self, i):
        top = math.pow(1.0+self.s, i) * GetXi(i)
        bottom = 0        
        for l in range(0,self.cellTypes+1):
            bottom += math.pow(1.0 +self.s, l)*GetXi(l)
        return top/bottom
        