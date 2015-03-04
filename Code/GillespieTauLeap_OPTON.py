# -*- coding: utf-8 -*-
"""
Created on Tue Feb 03 12:48:47 2015

@author: Connor
"""

import math
import random
import numpy
import gc

class BatchResult:
    def __init__(self):
        self.simCount = 0
        self.avgFixTime = 0.0
        self.avgFixProb = 0.0
        self.avgBadFrames = 0
        self.avgFrames = 0
        self.Reset()

    def Reset(self):
        self.simCount = 0
        self.avgFixTime = 0.0
        self.avgFixProb = 0.0
        self.avgBadFrames = 0
        self.avgFrames = 0
        

class Gillespie:
    def __init__(self):
        self.rateCallbacks = [] #Array of callbacks for all our rates
        self.eventCallbacks = []
        self.rateCallbackCount = 0
        self.rateCache = []
        
        self.eventCount = []
        
        self.params = 0
        self.lambd = 0.0
        self.simSteps = 0
        self.curTime = 0.0
        self.timeLimit = 10.0
        self.epsilon = 0.05        
        
        self.stopAtAppear = 1        
        
        self.history = 0
        
        self.preSim = 0
        self.postSim = 0
        
        self.tau = -1.0
        self.tauCache = []
        
        self.printProgress = 0 #Set to 1 to print batch progress. 2 prints batch and sim progress
        #Tau Debug Info
        self.RECORD_TAU_INFO = 0
        self.TAU_HIST = []
        self.BAD_FRAME_COUNT = 0
        self.BAD_FRAME_WARN = 0.05 #Fraction of frames that need to be bad to give a warning
        self.nextBatchProgressFrac = 0.1
        self.Reset()
        
    def Reset(self):
        self.simSteps = 0
        self.curTime = 0.0
        self.TAU_HIST = []
        self.BAD_FRAME_COUNT = 0
        self.nextProgressFrac = 0.1
        
    def AddCallback(self, rateFunc, stateChange):
        self.rateCallbacks.append(rateFunc)
        self.eventCallbacks.append(stateChange)
        self.rateCallbackCount += 1
        self.rateCache.append(0.0)
        self.eventCount.append(0)
        
    def Hook(self, param):
        self.params = param
        for i in range(0,len(self.params.n)):
            self.tauCache.append(0.0)
        param.Hook(self)
        
        print("GILLESPIE ADDED {0} CALLBACKS.".format(self.rateCallbackCount))
        
    def UnHook(self):
        self.rateCallbackCount = 0
        self.rateCallbacks = []
        self.eventCallbacks = []
        self.rateCache = []
            
    def Poission(self, mean):
        return numpy.random.poisson(mean)
        
    #Exponential parameter for frequency of events
    def GetLambda(self):
        return self.lambd

    #Get the axillary quanitity mu        
    def GetAuxMu(self, i):
        summation = 0.0        
        for j in range(0,self.rateCallbackCount):
            summation += self.eventCallbacks[j][i]*self.rateCache[j]
        #print("rate cache is: {0}".format(self.rateCache))
        #print("mu is given as: {0}".format(summation))
        if(summation == 0):
             summation = 0.001
             #print("WARNING: mu was found to be zero.")
             #print("rate cache is:".format(self.rateCache))
            
        return summation

    def GetAuxSigma(self, i):
        summation = 0.0        
        for j in range(0,self.rateCallbackCount):
            summation += (self.eventCallbacks[j][i] * self.eventCallbacks[j][i]) * self.rateCache[j]
        if(summation == 0):
            summation = 0.001
            #print("WARNING: sigma was found to be zero.")
        return summation
        
    #Returns an exponentially distributed number based on the lambda parameter
    def GetTimeStep(self):
        for i in range(0,len(self.params.n)):
            comp = max(self.epsilon*self.params.n[i], 1.0 )
            #print("mu is: {0} and sigma is: {1}".format(self.GetAuxMu(i),self.GetAuxSigma(i)) )            
            tau_element = min( comp / math.fabs(self.GetAuxMu(i)) , (comp * comp) / self.GetAuxSigma(i) )
            self.tauCache[i] = tau_element
            #print("tau is: {0}".format(self.tau))
            
        self.tau = min(self.tauCache)
        #self.tau = 0.01
        return self.tau
        
    #Chose and execute which event to carry out. Updates population counts
    #Uses weighted random number between 0 and 1
    def ChooseEvent(self):
        goodFrame  = 0
        while goodFrame == 0:        
            goodFrame = 1
            for i in range(0,self.rateCallbackCount):
                if self.rateCache[i] > 0.0001:
                    self.eventCount[i] = self.Poission(self.rateCache[i] * self.tau)
                else:
                    self.eventCount[i] = 0
                
                for pop in range(0, len(self.params.n)):
                    self.params.n[pop] += self.eventCallbacks[i][pop] * self.eventCount[i]
                    if self.params.n[pop] < 0:
                        goodFrame = 0
                        
            if goodFrame == 0:
                self.BAD_FRAME_COUNT += 1
                for i in range(0,self.rateCallbackCount):
                    for pop in range(0,len(self.params.n)):
                        self.params.n[pop] -= self.eventCallbacks[i][pop] * self.eventCount[i]
                        
                self.tau /= 2.0
  
    def UpdateRates(self):
        self.lambd = 0.0
        index = 0
        for i in range(0,self.params.typeCount):
            for j in range(0, self.params.typeCount):
                if i == j:
                    continue;
                self.rateCache[index] = self.params.GetTIJ(i,j)
                self.lambd += self.rateCache[index]
                
                index += 1
                
    def SetHistory(self, hist):
        self.history = hist
    
    def Simulate(self):
        gc.disable()
        self.params.Reset()
        self.Reset()
        
        if self.history != 0:
            self.history.RecordFrame(self.curTime, self.params)
        
        while self.curTime < self.timeLimit:
            if (self.printProgress >= 2 and float(self.curTime) / self.timeLimit >= self.nextSimProgressFrac):
                print(int(float(self.curTime) / self.timeLimit * 100.0)),
                self.nextSimProgressFrac += 0.25
                
            self.params.PreSim(self)
            if self.preSim != 0:
                self.preSim(self.curTime, self.params)
            self.UpdateRates()
            if self.lambd == 0 or (self.stopAtAppear == 1 and (self.params.n[self.params.typeCount - 1] >= 1)): #We have fixated at an absorbing state
                break
            self.tau = self.GetTimeStep() #How much time until the next event?
            self.ChooseEvent() #Choose what kind of event and update cell counts
            
            self.curTime += self.tau #Increment time
            if self.RECORD_TAU_INFO:
                self.TAU_HIST.append(self.tau)
            self.simSteps+= 1 #Increase event count
            
            self.params.PostSim(self)
            if self.postSim != 0:
                self.postSim(self.curTime, self.params)
            
            if self.history != 0:
                self.history.RecordFrame(self.curTime, self.params)
                
        if self.printProgress >= 2:
            print(' ')
                
        if float(self.BAD_FRAME_COUNT) / self.simSteps > self.BAD_FRAME_WARN:
            print("WARNING: {0}% Bad Frames. Consider lower epsilon!".format(int(float(self.BAD_FRAME_COUNT) / self.simSteps  * 100.0)))
        
        gc.enable()
                
    def SimulateBatch(self, simCount):
        res = BatchResult()
        res.simCount = simCount
        
        if simCount <= 0:
            return res        
        
        print("Batch %: 0"),
        for i in range(0, simCount):
            if (self.printProgress >= 1 and float(i) / simCount > self.nextBatchProgressFrac):
                print(int(float(i) / simCount * 100.0)),
                self.nextBatchProgressFrac += 0.1
            self.Simulate()
            if self.curTime < self.timeLimit: #we finished early
                res.avgFixProb += 1.0
                
            res.avgFixTime += self.curTime
            res.avgBadFrames += self.BAD_FRAME_COUNT
            res.avgFrames += self.simSteps

        if self.printProgress >= 1:
            self.nextBatchProgressFrac = 0.1
            print('100')
        
        res.avgFixProb = float(res.avgFixProb) / simCount
        res.avgFixTime = float(res.avgFixTime) / simCount
        res.avgBadFrames = float(res.avgBadFrames) / simCount
        res.avgFrames = float(res.avgFrames) / simCount
        return res
                