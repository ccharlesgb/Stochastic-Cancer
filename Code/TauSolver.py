# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 18:41:05 2015

@author: Connor
"""

import math

class Solver:
    def __init__(self, param):
        self.params = param
        self.CacheX0()
        #Newton-Raphson Params
        self.epsilon = 1e-14 #Min value f prime
        self.tolerance = 1e-9 # accuracy required
        self.maxIter = 200
        self.maxTau = 1e4
        self.tau_hist = []
        
    def ValidateJ(self, j):
        k = (self.params.typeCount - 1)
        if k == j:
            print("Getting out of bounds tau j = {0}".format(j))
        
    def CacheX0(self):
        self.params.Reset()
        self.X0_Cache = []
        for i in range(0, self.params.typeCount):
            self.X0_Cache.append(self.GetXj0(i))
        
    def GetXj0(self, j):
        if j == 0:
            return 1.0
        return max(self.params.GetU(j) * self.params.d * self.IntegralOfXj(j-1, 1.0), 1.0/self.params.N)
    
    def GetGammaTau(self, tau, j):
        tau = max(tau, 1e-99)
        s = self.params.r[1] - self.params.r[0]
        x_j_0 = self.X0_Cache[j]
        return math.sqrt(2.0/(s*tau) * math.log(1.0 / x_j_0))
    
    def IntegralOfXj(self,j, tau):
        if j < 0:
            return 0.0
        if j == 0:
            return tau * 1.0 #X_0(t) = 1.0
        
        s = self.params.r[1] - self.params.r[0]
        u = self.params.u[j]
        d = self.params.d
        
        
        gamma = self.GetGammaTau(tau, j)
        
        x_j_0 = self.X0_Cache[j]
        x_j_1 = self.X0_Cache[j-1]
        
        fac = 1.0/(s*gamma)
        selection = x_j_0 * (math.exp(s * gamma * tau)-1.0)
        mutation = u*d*x_j_1*(1.0/(s*gamma) - tau)        
        
        return fac * (selection + mutation)
        
    def GetXJ(self, t,j):
        if j < 0:
            return 0.0
        if j == 0:
            return 1.0           
        
        s = self.params.r[1] - self.params.r[0]
        u = self.params.GetU(j)
        d = self.params.d
        N = self.params.N   
        
        gamma = self.GetGammaTau(max(t,1e-10),j)        
        
        x_j_0 = self.X0_Cache[j]
        x_j_1 = self.X0_Cache[j-1]             
        
        fac = 1.0/(s*gamma)
        expFactor = 0.0*u*d*x_j_1 + x_j_0*s*gamma
        xj = fac * (expFactor * math.exp(s*gamma*t) - 0.0*u*d*x_j_1)
        return xj    
            
    #Original Model
    def GetTauOriginal(self):
        s = self.params.r[1] - self.params.r[0]
        numer = math.pow(math.log(s/(self.params.u[0]*self.params.d)),2.0)
        denom = 2.0 * s * math.log(self.params.N)
        return float(numer)/denom
    
    def GetWaitingTimeOriginal(self, k):
        if k == 21:
            print("off by one error! (probably did you mean celltypes - 1)")
        return k * self.GetTauOriginal()
        
    #Neglecting the transient phase
    def GetTauNeglect(self, j):
        s = self.params.r[1] - self.params.r[0]
        fac = (s/(self.params.GetU(j+1) * self.params.d))
        logs = math.log(1.0 + fac*math.sqrt(2.0*math.log(self.params.N)))
        top = math.pow(logs,2.0)
        bottom = 2.0 * s * math.log(self.params.N)
        return top/bottom
        
    def GetWaitingTimeNeglect(self, k):
        if k == 21:
            print("off by one error! (probably did you mean celltypes - 1)")
        j_i = -math.log(self.params.N)/math.log(self.params.u[0]*self.params.d)

        #Linear interpolate by taking the fractional part of the one we are
        #supposed to miss. ie if at 3.5 then skip the first 3 then skip half the fourth        
        frac_part = j_i - math.floor(j_i)
        first_whole = int(math.ceil(j_i))
        
        #Take the fractioanal part of the one we are missing
        frac_tau = frac_part * self.GetTauNeglect(int(math.floor(j_i)))
        
        total = 0.0
        for i in range(first_whole, k):
            total = total + self.GetTauNeglect(i)
        return total + frac_tau
        
        #TODO: Support variable mutation rate
        #return (k-j_i) * self.GetTauNeglect(0)
    
    
    #Modelling the transient phase with analytical
    def GetTauModel(self, j):
        self.ValidateJ(j)
        if j == 0:
            return 0.0 #Justified in maths
                
        s = self.params.r[1] - self.params.r[0]
        u = self.params.GetU(j+1)
        d = self.params.d
        N = self.params.N   
        x_j_0 = self.X0_Cache[j]
            
        twoLogx0 = math.sqrt(-2.0 * math.log(x_j_0))
        fac = 1.0 / (-2.0 * s * math.log(x_j_0))
        bracket = s/(N*u*d) * twoLogx0
        log = math.log(1.0 + bracket/x_j_0)
        
        return fac * math.pow(log,2.0)
    
    def GetWaitingTimeModel(self,k):
        if k == 21:
            print("off by one error! (probably did you mean celltypes - 1)")
        total = 0.0
        for i in range(0,k):
            total += self.GetTauModel(i)
        return total
    
    def GetWaitingTimeModelNew(self,k):
        if k == 21:
            print("off by one error! (probably did you mean celltypes - 1)")
        total = 0.0
        self.tau_hist = [0.0] * k
        for i in range(0,k):
            tau = self.GetTauModelNew(i)
            self.tau_hist.append(tau)
            total += tau
        return total
    
    def GetXJ_NEW(self,j,tau_j):
        s = self.params.r[1] - self.params.r[0]
        u = self.params.GetU(j+1)
        d = self.params.d
        N = self.params.N
        tau_j = max(tau_j, 1e-9)
        gamma = math.sqrt(2.0/(s*tau_j) * math.log(N))
        x_j_1 = 0.0
        if j == 1:
            x_j_1 = 1.0 #Assume all type 0 cells forever
        elif j > 1:
            x_j_1 = self.GetXJ_NEW(j-1,self.tau_hist[j-1])
        #x_j_1 = 0.0
        #print("TAU_J IS", tau_j)
        #print("Input into type {0} from type {1} = {2} TAU_-1 = {3}".format(j,j-1,x_j_1,self.tau_hist[j-1]))
        return 1.0 / (s * gamma) * ((u*d*x_j_1 + (s*gamma)/N)*math.exp(s*gamma*tau_j) - u*d*x_j_1)
    
    def GetTauModelNew(self,j):
        s = self.params.r[1] - self.params.r[0]
        up1 = self.params.GetU(j+1)
        d = self.params.d
        N = self.params.N
        if j == 0:
            self.tau_hist[j] = 1.0/(N * up1 * d)
            return self.tau_hist[j] #Neglect
            
        u = self.params.GetU(j)
        
        gamma = math.sqrt(2.0*math.log(N))
        x_j_1 = 0.0
        if j == 1:
            x_j_1 =1.0
        elif j > 1:
            x_j_1 = self.GetXJ_NEW(j-1,self.tau_hist[j-1])   
        #print("Solving tau_{0}, extra_term = {1}".format(j,x_j_1))
        #print("Term comparison {0} and {1}".format((s*gamma/N),u*d*x_j_1))
        fac = 1.0 / (2 * s * math.log(N))
        logx = ((s * gamma)**2)/(N * up1 * d)
        logx = logx / ((s*gamma/N) + u*d*x_j_1)
        result = fac * math.pow(math.log(1.0 + logx),2.0)
        self.tau_hist[j] = result
        #print("TAU_{0} = {1}".format(j, result))
        return result
            

        
    
    