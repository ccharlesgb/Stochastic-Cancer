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
        
    def RootFunc(self, j,tau):   #Find the roots of this function
        if j < 0:
            return 0.0
        s = self.params.r[1] - self.params.r[0]
        u = self.params.GetU(j)
        d = self.params.d
        N = self.params.N

        gamma = self.GetGammaTau(tau, j)
        x_j_0 = self.X0_Cache[j]
        x_j_1 = self.X0_Cache[j-1]
        
        selection = x_j_0 * (math.exp(s * gamma * tau)-1.0)
        mutation = 0.0
        if j > 0: #Avoid divide by zero this term is not valid for this case
            mutation = u*d*x_j_1*(1.0/(s*gamma) - tau)
            mutation2 = u*d*x_j_1*(1.0/(s*gamma))
            #print("DIFF = {0}".format(mutation - mutation2))
        offset = (s * gamma)/(N*u*d)
        return (selection + mutation) - offset
        
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
    
    #Newton Raphson to Solve tau
    def GetTau(self, j):
        self.ValidateJ(j)
        if j <= 2:
            return 0.0
        foundSol = False
        
        delta = 0.01
        x0 = 10.0
        s = self.params.r[1] - self.params.r[0]
        for i in range(0, self.maxIter):
            y = self.RootFunc(j, x0)
            y_prime = (self.RootFunc(j, x0 + delta) - y) / delta
            if abs(y_prime) < self.epsilon:
                break #Denominator too small
            x1 = min(x0 - y / y_prime, self.maxTau)
            
            if x1 == 0  or abs(x1 - x0)/abs(x1) < self.tolerance:
                foundSol = True
                break
            x0 = x1
        if foundSol == True:
            return x0
        else:
            return -1.0
            
    def GetWaitingTime(self, k):
        total = 0.0
        for i in range(0, k):
            total += self.GetTau(i)
        return total
            
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
        x_j_1 = 0.0
        if j > 0:
            x_j_1 = self.X0_Cache[j-1]
            
            
            
        twoLogx0 = math.sqrt(-2.0 * math.log(x_j_0))
        fac = 1.0 / (-2.0 * s * math.log(x_j_0))
        #bracket = s/(N*u*d) * twoLogx0 - (0.0* u * d * x_j_1)/(s * twoLogx0)
        bracket = s/(N*u*d) * twoLogx0
        #bracket = max(bracket, 0.0)    
        
        log = math.log(1.0 + bracket/x_j_0)
        
        return fac * math.pow(log,2.0)
    
    def GetWaitingTimeModel(self,k):
        if k == 21:
            print("off by one error! (probably did you mean celltypes - 1)")
        total = 0.0
        for i in range(0,k):
            total += self.GetTauModel(i)
        return total
            

        
    
    