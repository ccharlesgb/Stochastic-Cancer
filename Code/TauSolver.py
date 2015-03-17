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
        return max(self.params.u[0] * self.params.d * self.Xj_Integral(j-1, 1.0), 1.0/self.params.popSize)

    def Xj_Integral(self, j, tau):
        s = self.params.r[1] - self.params.r[0]
        u = self.params.u[0]
        d = self.params.d
        N = self.params.popSize        
        u = 0.0
        if j < 0:
            return 0.0
        if j == 0:
            return tau
        x_j_0 = self.X0_Cache[j]
        x_j_1 = self.X0_Cache[j-1]
        
        gamma = math.sqrt(2.0 / (s * tau) * math.log(1.0 / (x_j_1)))
        #a = s * gamma - u * d
        b = u * d * x_j_1
        c = x_j_0
        
        a = s * gamma
        
        #top = (a*c + b)*math.exp(a * tau) - a * b * tau - a *c - b
        top = (a*c + b)*math.exp(a * tau) - a * b * tau - a * c
        top = top / (a*a)
        top = top
        
        return top 

    def GetTau(self, j):
        foundSol = False
        
        delta = 0.01
        x0 = 10.0
        c = 1.0 / (self.params.popSize * self.params.u[0] * self.params.d)
        for i in range(0, self.maxIter):
            y = self.Xj_Integral(j, x0) - c
            y_prime = ((self.Xj_Integral(j, x0 + delta) - c) - y) / delta
            #print("x0 = {0} y' = {1}".format(x0, y_prime))
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
            
    def GetXJ(self, t,j):
        s = self.params.r[1] - self.params.r[0]
        u = self.params.u[0]
        d = self.params.d
        N = self.params.popSize    
        
        x_j_0 = 1.0 / N
        x_j_1 = 1.0 / N     
        if j == 0:
            x_j_0 = 1.0 - 1.0/N
            x_j_1 = 0.0
        if j == 1:
            x_j_0 = 1e4 / N
            x_j_1 = 1.0 - 1.0/N
        if j == 2:
            x_j_1 = 1e4 / N 
        t = float(t)
        tau_epsilon = 1e-10
        gamma = math.sqrt((2.0 / (s * (t+tau_epsilon))) * math.log(1.0 / (x_j_1)))
        #gamma = math.sqrt(2 * math.log(1.0 / (x_j_0)))
        a = s * gamma
        b = u * d * x_j_1
        c = x_j_0
        a = s * gamma
        
        res = (b + c*a)*math.exp(a*t) - b
        return 1.0/a * res
        
    
    