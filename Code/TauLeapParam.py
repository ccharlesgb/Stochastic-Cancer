# -*- coding: utf-8 -*-
"""
Created on Mon Mar 09 14:41:29 2015

@author: Connor
"""
import math
import scipy.misc as scimisc
import hashlib

class Hist:
    def __init__(self, typeCount):
        self.typeCount = typeCount
        self.ClearFrames()
        self.ResetTotal()
        self.SDP = 1

    def ResetTotal(self):
        self.tTotal = []
        self.histArrayTotal = dict()
        for i in range(0, self.typeCount):
            self.histArrayTotal[i] = []
       
    def ClearFrames(self):
        self.histArray = dict()
        self.tHist = []
        self.yearHist = []
        self.thetajHist = dict()
        self.avgJHist = []
        self.avgSJHist = []
        for i in range(0, self.typeCount):
            self.histArray[i] = []
            self.thetajHist[i] = []

    def RecordFrame(self, sim):
        self.tHist.append(sim.curTime)
        self.yearHist.append(float(sim.curTime) / 365.0)
        totalJ = 0.0
        totalSJ = 0.0            
        if len(self.tHist) >= len(self.tTotal): #Append average hist if neccessary
            self.tTotal.append(sim.curTime)
        for i in range(0, self.typeCount):
            self.histArray[i].append(sim.n[i])
            
            log_ni = 0.0
            if sim.n[i] > 0:
                log_ni = math.log10(float(sim.n[i])) / self.SPD
            if len(self.histArray[i]) >= len(self.histArrayTotal[i]): #Append average hist if neccessary
                self.histArrayTotal[i].append(log_ni)
            else:
                self.histArrayTotal[i][sim.curTime] += log_ni
            
            if sim.prob_vector != 0:
                self.thetajHist[i].append(sim.prob_vector[i])
            totalJ += i * float(sim.n[i]/sim.params.N)
            totalSJ += (sim.params.r[i] - 1.0 )* float(sim.n[i]/sim.params.N)
        self.avgJHist.append(totalJ)
        self.avgSJHist.append(totalSJ)
         
    def GetDictionary(self):
        runDict = dict()
        return runDict

#Define all the parameters for our stem cell model
class Params:
    def __init__(self, typeCount):
        self.SetTypeCount(typeCount)
        self.d = 100
        self.USE_D = True
        self.avgFit = 0.0
        self.uNotConst = 0
        self.Reset()
        
        self.hash = hashlib.md5()
        
    def SetTypeCount(self, typeCount):
        self.typeCount = typeCount
        self.n0 = []
        self.r = []
        self.u = []
        self.thetaJCache = []

        for i in range(0, typeCount):
            self.n0.append(0)
            self.r.append(1.0)
            self.thetaJCache.append(0.0)
        for i in range(0, typeCount-1):
            self.u.append(0.1)
        self.n0[0] = 1e2
        
    #Returns a hash of all the parameters so you can tell if theyre the same
    def GetHash(self):
        for i in range(0,self.typeCount):
            self.hash.update(str(self.typeCount))
            self.hash.update(str(self.n0[i]))
            self.hash.update(str(self.r[i]))
            self.hash.update(str(self.u[i]))
        self.hash.update(str(self.d))
        self.hash.update(str(self.uNotConst))
        self.hash.update(str(self.USE_D))
        return self.hash.hexdigest()
        
    #Returns a nice file string to summarise parameters
    def GetFileString(self):
        N_format = str(self.N)
        if self.N >= 1e3:
            power = math.log(self.N, 10)
            val = round(self.N / (10.0**power),1) #Round to 1 DP
            N_format = "{0}e{1}".format(val,power)
        myFileStr = "K_{0}_N_{1}_d_{2}_u_{3}_s_{4}".format(self.typeCount, N_format, self.d, self.u[0], self.r[1] - self.r[0])
        return myFileStr
        
    def Reset(self):
        self.N = 0        
        for pop in range(0,self.typeCount):
            self.N += self.n0[pop]
        self.CacheCombinations()
            
    def EventTIJ(self, i,j):
        arr = []
        for x in range(0,self.typeCount):
            arr.append(0)
            if x == i:
                arr[x] = -1
            elif x == j:
                arr[x] = 1
        return arr

    #Todo calculate this once
    def GetAvgFit(self, n):
        tot = 0.0
        for i in range(0, self.typeCount):
            tot += self.r[i] * n[i]
        return tot
    
    def Hook(self, sim):
        for i in range(0,self.typeCount):
            for j in range(0, self.typeCount):
                if i != j:
                    sim.AddStateChange(i, j, self.EventTIJ(i,j))
                    #print("Adding rate {0},{1}".format(i,j))
    
    #Reaction probability for cell from 1->0
    def GetTIJ(self, i, j, n):
        return self.thetaJCache[j] * n[i]# * math.sqrt(2.0)
        
    def SetUAll(self,u):
        for i in range(1,self.typeCount):
            self.SetU(i,u)

    def SetU(self, j, u):
        if j <= 0:
            print("Invalid u_{0} too small".format(j))
        if j > self.typeCount - 1:
            print("Invalid u_{0} index too large".format(j))
        self.u[j-1] = u
        
    def GetU(self, j):
        if j <= 0:
            print("Invalid j < 1")
        if j > self.typeCount - 1:
            print("Invalid j > k")
        
        u = self.u[j-1]
        return u
            

    #Because T_i_j is just n_i * theta_j can just do this once
    #Much faster, not sure if this is identical to the other theta j?
    #Only used for tau leap not wright fisher need to merge them
    def CacheThetaJ(self, n):
        for j in range(0, self.typeCount):
            if n[j] == 0 and (j > 0 and n[j-1] == 0): #If we dont have any cells we cant divide, if we dont have any previous cells we cant get mutated to
                self.thetaJCache[j] = 0.0
                continue
            top = 0.0
            #u_j = self.u[j] * (self.d - j)
            #u_jm1 = self.u[j-1] * (self.d - (j-1))
            u_j = 0.0
            if j < self.typeCount - 1:
                u_j = self.GetU(j+1) #Mutation rate from type j to j+1
            u_jm1 = 0.0
            
            if self.USE_D == True: #Use the susceptible loci model?
                u_j  = u_j * (self.d - j)
                u_jm1 = u_jm1 * (self.d - j + 1)
                
            if j > 0: #No such thing as u_0
                u_jm1 = self.GetU(j) #Mutation rate from type j-1 to j
            
            if j == 0:
                top = (self.r[j] * (1.0 - u_j) * n[j])
            elif j == self.typeCount - 1:
                top = (self.r[j]*n[j] + self.r[j-1] * u_jm1 * n[j-1])
            else:
                top = (self.r[j]*(1.0 - u_j)*n[j] + self.r[j-1]*u_jm1*n[j-1])
            rate = top / self.avgFit
            self.thetaJCache[j] = rate
    
    def SetCompoundFitness(self,s):
        for i in range(0,self.typeCount):
            self.r[i] = math.pow(1.0 + s, i)    
    
    def CacheCombinations(self):
        self.combinations = [[]]
        for i in range(0,self.typeCount):
            self.combinations.append([])
            for j in range(0,self.typeCount):
                comb = scimisc.comb(self.d-i, j-i)
                self.combinations[i].append(comb)        
    
    #Mutation?
    #u_0 = is invalid surely?
    #u_j = mutation rate from type j-1 to type j
    #There is a u_20 if celltypes = 21
    #But this is indexed as u[19]
    def GetThetaj(self,j, n):
        #if self.uNotConst == 0:
        summation = 0.0
        avgFit = self.avgFit
        if self.USE_D == False:
            for i in range(0, j+1):
                summation += math.pow(self.GetU(i+1), j-i)*math.pow(1.0-self.GetU(i+1), self.typeCount - 1 - j)*(self.r[i] * n[i])/avgFit 
        else:
            for i in range(0, j+1):
                summation += self.combinations[i][j]*math.pow(self.GetU(i+1), j-i)*math.pow(1.0-self.GetU(i+1), self.d-j)*(self.r[i] * n[i])/avgFit
        return summation

    def PreSim(self, gillespie):
        self.avgFit = self.GetAvgFit(gillespie.n)
        self.CacheThetaJ(gillespie.n)
    
    def PostSim(self, gillespie):
        return
