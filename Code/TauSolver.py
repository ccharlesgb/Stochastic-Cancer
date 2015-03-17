# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 14:59:13 2015

@author: Connor
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 18:41:05 2015

@author: Connor
"""

import math
import matplotlib.pyplot as plt
import wright_fisher

class Solver:
    def __init__(self, param):
        self.params = param
        self.CacheX0()
        #Newton-Raphson Params
        self.epsilon = 1e-14 #Min value f prime
        self.tolerance = 1e-5 # accuracy required
        self.maxIter = 200
        self.maxTau = 1e4
        
    def CacheX0(self):
        self.X0_Cache = []
        for i in range(0, self.params.cellTypes + 1):
            self.X0_Cache.append(self.GetXj0(i))
        
    def GetXj0(self, j):
        return max(self.params.u[0] * self.params.d * self.IntegralOfXj(j-1, 1.0), 1.0/self.params.popSize)
    
    def GetGammaTau(self, tau):
        tau = max(tau, 1e-99)
        s = self.params.r[1] - self.params.r[0]
        N = self.params.popSize         
        return math.sqrt(2.0/(s*tau) * math.log(N)) 
    
    def IntegralOfXj(self,j, tau):
        s = self.params.r[1] - self.params.r[0]
        u = self.params.u[0]
        d = self.params.d
        
        if j < 0:
            return 0.0
        if j == 0:
            return tau        

        gamma = self.GetGammaTau(tau)      
        
        x_j_0 = self.X0_Cache[j]
        x_j_1 = self.X0_Cache[j-1]        
        
        fac = 1.0/(s*gamma)
        selection = x_j_0 * (math.exp(s * gamma * tau)-1.0)
        mutation = u*d*x_j_1*(1.0/(s*gamma) - tau)        
        
        return fac * (selection + mutation)

    #Newton Raphson to Solve tau
    def GetTau(self, j):
        foundSol = False
        
        delta = 0.01
        x0 = 10.0
        s = self.params.r[1] - self.params.r[0]
        for i in range(0, self.maxIter):
            c = (1.0) / (self.params.u[0] * self.params.d * self.params.popSize)
            print(c)
            y = self.IntegralOfXj(j, x0) - c
            y_prime = ((self.IntegralOfXj(j, x0 + delta) - c) - y) / delta
            print("j ={0} x0 = {1} y' = {2}".format(j, x0, y_prime))
            if abs(y_prime) < self.epsilon:
                break #Denominator too small
            x1 = min(x0 - y / y_prime, self.maxTau)
            
            if abs(x1 - x0)/abs(x1) < self.tolerance:
                foundSol = True
                break
            x0 = x1
        if foundSol == True:
            return x0
        else:
            return -1.0
            
    def GetWaitingTime(self, k):
        total = 0.0
        for i in range(0, k + 1):
            total += self.GetTau(i)
        return total
        
    def GetWaitingTimeOriginal(self, k):
        s = self.params.r[1] - 1.00
        numer = k * math.pow(math.log(s/(self.params.u[0]*self.params.d)),2.0)
        denom = 2.0 * s * math.log(self.params.popSize)
        return float(numer)/denom
    
    def GetWaitingTimeNeglect(self, k):
        j_i = -math.log(self.params.popSize)/math.log(self.params.u[0]*self.params.d)
        print("J_I", j_i)
        return (k-j_i) * self.GetTauNeglect()
            
    def GetXJ(self, t,j):
        if j < 0:
            return 0.0
        if j == 0:
            return 1.0           
        
        s = self.params.r[1] - self.params.r[0]
        u = self.params.u[0]
        d = self.params.d
        N = self.params.popSize    
        
        gamma = self.GetGammaTau(max(t,1e-10))        
        
        x_j_0 = self.X0_Cache[j]
        x_j_1 = self.X0_Cache[j-1]             
        
        fac = 1.0/(s*gamma)
        expFactor = u*d*x_j_1 + x_j_0*s*gamma
        xj = fac * (expFactor * math.exp(s*gamma*t) - u*d*x_j_1)
        return xj
        
        
    def GetTauNeglect(self):
        s = self.params.r[1] - self.params.r[0]
        logs = math.log(1.0 + s / (self.params.u[0] * self.params.d)*math.sqrt(2.0*math.log(self.params.popSize)))
        top = math.pow(logs,2.0)
        bottom = 2.0 * s * math.log(self.params.popSize)
        return top/bottom
        
    
    