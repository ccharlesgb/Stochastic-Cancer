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
        self.timeLimit = 100.0
        
        self.history = 0
        
        self.preSim = 0
        self.postSim = 0
        
    def AddCallback(self, rateFunc, eventFunc):
        self.rateCallbacks.append(rateFunc)
        self.eventCallbacks.append(eventFunc)
        self.rateCallbackCount += 1
        self.rateCache.append(0.0)
        
    def Hook(self, param):
        self.params = param
        param.Hook(self)
        
    def UnHook(self):
        self.rateCallbackCount = 0
        self.rateCallbacks = []
        self.eventCallbacks = []
        self.rateCache = []
        
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
        for i in range(0,self.rateCallbackCount):
            threshold += self.rateCache[i]
            if (rand < threshold):
                self.eventCallbacks[i]()
                return
            
    def UpdateRates(self):
        self.lambd = 0.0
        for i in range(0,self.rateCallbackCount):
            self.rateCache[i] = self.rateCallbacks[i]()
            self.lambd += self.rateCache[i]
        #print("SELF.LAMBDA = ", self.lambd)
            
    def Reset(self):
        self.simSteps = 0
        self.curTime = 0.0       
    
    def SetHistory(self, hist):
        self.history = hist
    
    def Simulate(self):
        self.params.Reset()
        self.Reset()
        
        if self.history != 0:
            self.history.RecordFrame(self.curTime, self.params)
        
        while self.curTime < self.timeLimit:   
            self.params.PreSim(self)
            if self.preSim != 0:
                self.preSim(self.curTime, self.params)
                
            self.UpdateRates()
            if self.lambd == 0: #We have fixated at an absorbing state
                print("END")
                return
                
            timestep = self.GetTimeStep() #How much time until the next event?
            self.ChooseEvent() #Choose what kind of event and update cell counts
            
            self.curTime += timestep #Increment time
            self.simSteps+= 1 #Increase event count
            
            self.params.PostSim(self)
            if self.postSim != 0:
                self.postSim(self.curTime, self.params)
            
            if self.history != 0:
                self.history.RecordFrame(self.curTime, self.params)
                