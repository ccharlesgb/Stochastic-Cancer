# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 03 12:01:50 2015

@author: Connor
"""

import random

class NHist:
    def __init__(self, typeCount):
        self.typeCount = typeCount
        self.ClearFrames()
       
        
    def ClearFrames(self):
        self.tHist = []
        self.histArray = dict()
        for i in range(0, self.typeCount):
            self.histArray[i] = []
            
        
        
    def RecordFrame(self, time, param):
        if random.random() < 1.1:
            self.tHist.append(time)
            for i in range(0, self.typeCount):
                self.histArray[i].append(param.ni[i])
         
    def GetDictionary(self):
        runDict = dict()

        
        return runDict

#Define all the parameters for our stem cell model
class NParam:
    def __init__(self, typeCount):
        #Define Population Counts
        self.ni = []
        self.ini = []
        self.ri = []
        self.ui = []
        self.typeCount = typeCount
        for i in range(0,self.typeCount):
            self.ni.append(0)
            self.ini.append(0)
            self.ri.append(float(i) * 0.02)
            self.ui.append(0.001)
        self.Reset()
            
    def Reset(self):
        for i in range(0,self.typeCount):
            self.ni[i] = self.ini[i]
        
    
    def PreSim(self, gillespie):
        return
        
    def PostSim(self, gillespie):
        return     
        