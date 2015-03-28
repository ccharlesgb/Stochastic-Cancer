# -*- coding: utf-8 -*-
"""
Created on Tue Feb 03 12:48:47 2015

@author: Connor
"""

import math
import random

class Gillespie:
    def __init__(self):
        self.rateCallbacks = [] #Array of callbacks for all our rates
        self.eventCallbacks = []
        self.rateCallbackCount = 0
        self.rateCache = []
        
        self.params = 0
        self.lambd = 0.0
        self.simSteps = 0
        self.curTime = 0.0
        self.timeStep = 0.0
        self.timeLimit = 100.0
        
        self.history = 0
        
        self.preSim = 0
        self.postSim = 0
        
        self.changeEventsTot = 0.0
        self.mutEventsTot = 0.0
        
    #Exponential parameter for frequency of events
    def GetLambda(self):
        return self.lambd

    #Returns an exponentially distributed number based on the lambda parameter
    def GetTimeStep(self):
        return 1.0/self.GetLambda() * math.log(1.0/random.random())
        
    #Chose and execute which event to carry out. Updates population counts
    #Uses weighted random number between 0 and 1
    def ChooseEvent(self):
        rand = random.random() * self.lambd
        
        threshold = 0.0
        cacheDex = 0
        for init in range(0, self.params.typeCount):
            for final in range(0, self.params.typeCount):
                if init == final:
                    continue
                #print("INIT", init, "FINAL", final, "RATE", self.rateCache[cacheDex])
                threshold += self.rateCache[cacheDex]
                cacheDex += 1
                if (rand < threshold):
                    self.params.ni[init] -= 1
                    self.params.ni[final] += 1
                    return
        for i in range(0,self.params.typeCount - 1):
            #print("MUT RATE", i, self.rateCache[cacheDex])
            threshold += self.rateCache[cacheDex]
            cacheDex += 1
            if (rand < threshold):
                self.params.ni[i] -= 1
                self.params.ni[i+1] += 1
                return
            
    def UpdateRates(self):
        self.lambd = 0.0
        cacheDex = 0
        self.rateCache = []
        self.changeEventsTot = 0.0
        for init in range(0, self.params.typeCount):
            for final in range(0,self.params.typeCount):
                if init == final:
                    continue
                self.rateCache.append(self.GetRate(init, final))
                self.lambd += self.rateCache[cacheDex]
                self.changeEventsTot += self.rateCache[cacheDex]
                cacheDex += 1
                
        self.mutEventsTot = 0.0
        for mut in range(0, self.params.typeCount-1):
            self.rateCache.append(self.params.ui[mut]*self.params.ni[mut])
            self.mutEventsTot += self.rateCache[cacheDex]
            self.lambd += self.rateCache[cacheDex]
            cacheDex += 1
            
        
    def Reset(self):
        self.simSteps = 0
        self.curTime = 0.0      
        self.timeStep = 0.0
    
    def SetHistory(self, hist):
        self.history = hist
    
    def GetRate(self, init, final):
        averageFit = 0.0
        
        for i in range(0, self.params.typeCount):
            averageFit += (1.0 + self.params.ri[i]) * self.params.ni[i]
    
        rate = self.params.ni[init] * (1.0 + self.params.ri[final])*self.params.ni[final]
        rate = float(rate) / float(averageFit)
        return rate
    
    def Simulate(self):
        self.params.Reset()
        self.Reset()
        
        if self.history != 0:
            self.history.RecordFrame(self.curTime, self.params)
      
        while self.curTime < self.timeLimit:   
            #self.params.PreSim(self)
            if self.preSim != 0:
                self.preSim(self.curTime, self.params)
                
            self.UpdateRates()
            if self.lambd == 0: #We have fixated at an absorbing state
                return
                
            self.timeStep = self.GetTimeStep() #How much time until the next event?
            self.ChooseEvent() #Choose what kind of event and update cell counts
            
            self.curTime += self.timeStep #Increment time
            self.simSteps+= 1 #Increase event count
            
            self.params.PostSim(self)
            if self.postSim != 0:
                self.postSim(self.curTime, self.params)
            
            if self.history != 0:
                self.history.RecordFrame(self.curTime, self.params)
                