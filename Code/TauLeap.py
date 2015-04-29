# -*- coding: utf-8 -*-
"""
Created on Mon Mar 09 14:07:50 2015

@author: Connor
"""

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
        

class Sim:
    def __init__(self, typeCount):
        self.typeCount = typeCount        
        self.typeRange = range(0,typeCount)        
        
        self.eventIDs = []          #i->j tuples
        self.stateChange = []        #Array of state change vectors
        self.rateCount = 0        #Number of reactions that are possible
        self.rateCache = []        #A cache of the rates for the current state.
        self.eventCount = []    #How many times an event will fire        
        
        self.params = 0
        
        self.n = [] #Our state vector
        for i in self.typeRange:
            self.n.append(0)
        
        self.critLambda = 0.0
        self.lambd = 0.0
        self.simSteps = 0
        self.curTime = 0.0
        self.timeLimit = 10.0
        
        self.epsilon = 0.05 #Error parameter
        self.n_c = 10 #How close a reaction is to depleting its resource before its considered critical
        self.EnableTauLeap = True #If tau leaping is disabled all reactions are critical (Gillespie)
        self.w = 10 #If our dynamic tau < w / self.lambd then abondon tau leaping (its bogged down in critical reactions)
        
        
        self.stopAtAppear = 1        
        
        self.history = 0
        
        self.preSim = 0
        self.postSim = 0
        
        self.tau_prime = -1.0 #Candidate tau from dynami selection
        self.tau_prime2 = -1.0 #Candidate tau from critical set of events
        
        self.criticalCount = 0
        self.criticalRateIndices = []
        self.criticalTypeIndices = [False] * typeCount
        
        self.tau = -1.0
        self.tauCache = []
        
        self.prob_vector = 0 #Compatibility
        
        self.MAX_POISSION = 2.1e9 #Max poission number        
        
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
        self.critFrames = 0
        self.BAD_FRAME_COUNT = 0
        self.nextProgressFrac = 0.1
        
        self.rateRange = range(0,self.rateCount)        
        
        if self.history != 0:
            self.history.ClearFrames()
        
        if self.params != 0:
            for i in self.typeRange:
                self.n[i] = self.params.n0[i]
    
    def AddStateChange(self, i, j, stateChange):
        self.eventIDs.append([i,j])
        self.stateChange.append(stateChange)
        self.rateCount += 1
        self.rateCache.append(0.0)
        self.eventCount.append(0)
        self.criticalRateIndices.append(False) #False is non critial true is critical
    
    def Poission(self, mean):
        return numpy.random.poisson(mean)
  
    def UpdateCriticalSet(self):
        #Loop through all our rates
        self.criticalCount = 0
        self.criticalTypeCount = 0
        self.criticalTypeIndices = [False] * self.typeCount;
        for i in self.rateRange:
            self.criticalRateIndices[i] = False #Assume it isnt critical
            if self.rateCache[i] > 0.0: #Can this rate even happen?
                i_j = self.eventIDs[i] #Get the i->j this corresponds to
                if self.EnableTauLeap == False or self.n[i_j[0]] < self.n_c: #Is the cell type that this event depletes (ie 'i') less than the critical number?
                    self.criticalRateIndices[i] = True
                    self.criticalCount += 1
                    self.criticalTypeIndices[i_j[0]] = True
        #print("Critical indices", self.criticalIndices, self.n)
                
            
    def UpdateRates(self):
        self.lambd = 0.0
        for i in self.rateRange:
            i_j = self.eventIDs[i]
            if self.n[i_j[0]] == 0:
                self.rateCache[i] = 0.0
            else:
                self.rateCache[i] = self.params.GetTIJ(i_j[0],i_j[1], self.n)
            if self.rateCache[i] < 0.0:
                print("RATE IS NEGATIVE i,j = {0} {1}".format(*i_j))
                self.curTime = self.timeLimit
            #print("Getting rate ", i_j, self.rateCache[i])
            self.lambd += self.rateCache[i]

    #Get the mean change in the population type i
    def GetMeanChange(self, i):
        summation = 0.0        
        for j in self.rateRange:
            if self.criticalRateIndices[j] == False:
                summation += self.stateChange[j][i]*self.rateCache[j]
        if(summation == 0.0):
             summation = 1e-99
        return summation
        
    #Get the standard deviation of the population type i
    def GetSDChange(self, i):
        summation = 0.0
        for j in self.rateRange:
            if self.criticalRateIndices[j] == False:
                summation += (self.stateChange[j][i] * self.stateChange[j][i]) * self.rateCache[j]
        if(summation == 0.0):
            summation = 1e-99
        return summation

    def UpdateTauPrime(self):
        if self.criticalCount == self.rateCount: #There are no non-critical reactions
            self.tau_prime = 1e99
            return
        
        self.tau_prime = 1e99
        tau_element = 1e99
        for i in self.typeRange:
            if self.criticalTypeIndices[i] == False: #Is this indice critical?
                toleratedChange = max(self.epsilon*self.n[i], 1.0 )
                    
                meanChange  = self.GetMeanChange(i)
                sdChange = self.GetSDChange(i)
                tau_element = min( toleratedChange / math.fabs(meanChange) , (toleratedChange * toleratedChange) / sdChange )
                #print("Type {0} Tolerated {1} Mean {2} SD {3} Element {4}".format(i, toleratedChange, meanChange, sdChange, tau_element))
                self.tau_prime = min(self.tau_prime, tau_element)
        
    def UpdateTauPrime2(self):
        if self.criticalCount == 0: #There are no critical reactions
            self.tau_prime2 = 1e99
            return
        
        self.critLambda = 0.0
        for i in self.rateRange:
            if self.criticalRateIndices[i] == True:
                self.critLambda += self.rateCache[i]
        self.tau_prime2 = 1.0/self.critLambda * math.log(1.0/random.random())
        #return random.expovariate(self.GetLambda())
        
    def ExecuteCrit(self):
        r_1 = random.random() * self.critLambda #Number between 0 and 1
        threshold = 0.0
        for i in self.rateRange:
            if self.criticalRateIndices[i] == True:
                threshold += self.rateCache[i]
                #print("r_1 = {0} thresh = {1}".format(r_1, threshold))
                if r_1 < threshold:
                    for pop in self.typeRange:
                        self.n[pop] += self.stateChange[i][pop]
                    #print("Executed rate", i)
                    return
    
    def ExecuteNonCrit(self):
        goodFrame = 0
        while goodFrame == 0:
            goodFrame = 1
            for i in self.rateRange:
                if self.criticalRateIndices[i] == False:
                    if self.rateCache[i] > 1e-10:
                        mean = self.rateCache[i] * self.tau
                        
                        if mean >= self.MAX_POISSION:
                            print("Dividing Poission Mean from {0}".format(mean))
                            gen_count = 10
                            p_total = 0.0
                            for i in range(0, gen_count):
                                p_total += self.Poission(float(mean)/gen_count)
                                
                            self.eventCount[i] = p_total
                        else:
                            self.eventCount[i] = self.Poission(mean)
                    else:
                        self.eventCount[i] = 0
                        continue
                    for pop in self.typeRange:
                        self.n[pop] += self.stateChange[i][pop] * self.eventCount[i]
                        
            for pop in self.typeRange:
                if self.n[pop] < 0:
                    #print("BAD FRAME TRUE BECAUSE n[{0}] = {1}".format(pop, self.n[pop]))
                    goodFrame = 0
                                
            if goodFrame == 0:
                print("BAD FRAME", self.n)
                self.BAD_FRAME_COUNT += 1
                for i in self.rateRange:
                    if self.criticalRateIndices[i] == False:
                        for pop in self.typeRange:
                            self.n[pop] -= self.stateChange[i][pop] * self.eventCount[i]
                        
                self.tau /= 2.0
    
    def Simulate(self):
        gc.disable()
        self.params.Reset()
        self.Reset()
        
        if self.history != 0:
            self.history.RecordFrame(self)
        
        while self.curTime < self.timeLimit:
            #Print Progress
            #print("N IS ", self.n)
            if (self.printProgress >= 2 and float(self.curTime) / self.timeLimit >= self.nextSimProgressFrac):
                print(int(float(self.curTime) / self.timeLimit * 100.0)),
                self.nextSimProgressFrac += 0.25
            #Do Pre sim callbacks
            self.params.PreSim(self)
            if self.preSim != 0:
                self.preSim(self.curTime, self.params)
            #Start the actual sim
            #Update our rates for each reaction
            self.UpdateRates()
            self.UpdateCriticalSet() #Update the reactions that are critical
            
            #Quick check for fixation or appearance of the mutant we want
            if self.lambd == 0 or (self.stopAtAppear == 1 and (self.n[self.typeCount - 1] >= 1)):
                print("Stopped because appeared", self.stopAtAppear)
                break
            
            self.UpdateTauPrime()
            self.UpdateTauPrime2()
            
            self.tau = self.tau_prime #Assume no critical reactions
            #print(self.tau_prime, self.tau_prime2)
            if self.tau_prime >= self.tau_prime2: #Chose a critical reaction to occur
                #print("Crit Reaction")
                self.critFrames += 1
                self.tau = self.tau_prime2 #Use the critical waiting time instead
                self.ExecuteCrit()
            self.ExecuteNonCrit() #Always execute non critical reactions
            
            self.curTime += self.tau #Increment time
            if self.RECORD_TAU_INFO:
                self.TAU_HIST.append(self.tau)
            self.simSteps+= 1 #Increase event count
            
            self.params.PostSim(self)
            if self.postSim != 0:
                self.postSim(self.curTime, self.params)
            
            if self.history != 0:
                self.history.RecordFrame(self)
                
        if self.printProgress >= 2:
            print(' ')
        #print("Crit Frames: {0}%".format(round(float(self.critFrames) / self.simSteps * 100.0)))  
        if self.simSteps != 0 and float(self.BAD_FRAME_COUNT) / self.simSteps > self.BAD_FRAME_WARN:
            print("WARNING: {0}% Bad Frames. Consider lower epsilon!".format(int(float(self.BAD_FRAME_COUNT) / self.simSteps  * 100.0)))
        
        gc.enable()
                
    def SimulateBatch(self, simCount):
        res = BatchResult()
        res.simCount = simCount
        
        if simCount <= 0:
            res.avgFixTime = 0.0
            return res        
        
        if self.printProgress >= 1:
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
                