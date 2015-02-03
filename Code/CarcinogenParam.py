# -*- coding: utf-8 -*-
"""
Created on Tue Feb 03 13:10:36 2015

@author: Connor
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Feb 03 12:01:50 2015

@author: Connor
"""

#Define all the parameters for our stem cell model
class CarcinogenParam:
    def __init__(self):
        self.n0 = 0
        self.n1 = 0
        self.n2 = 0
        
        self.in0 = 10
        self.in1 = 0
        self.in2 = 0
        
        self.r0 = 1.0
        self.r1 = 1.0
        self.r2 = 1.0
        
        self.u1 = 0.1
        self.u2 = 0.1
        
    def AvgFitness(self):
        return (self.r0*self.n0 + self.r1*self.n1 + self.r2 * self.n2)
        
    def Reset(self):
        self.n0 = self.in0
        self.n1 = self.in1
        self.n2 = self.in2
        
        self.N = self.in0 + self.in1 + self.in2
    
    def Hook(self, gillespie):
        gillespie.AddCallback(self.GetT10, self.EventT10)
        gillespie.AddCallback(self.GetT01, self.EventT01)
        gillespie.AddCallback(self.GetT20, self.EventT20)
        gillespie.AddCallback(self.GetT02, self.EventT02)
        gillespie.AddCallback(self.GetT12, self.EventT12)
        gillespie.AddCallback(self.GetT21, self.EventT21)
    
    #Reaction probability for cell from 1->0
    def GetT10(self):
        top = (1.0 - self.u1) * self.r0 * float(self.n0) * self.n1
        return top / self.AvgFitness()
    
    def GetT20(self):
        top = (1.0 - self.u1) * self.r0 * float(self.n0) * self.n2
        return top / self.AvgFitness()
    
    def GetT01(self):
        top = (self.u1 * self.r0 * self.n0 + (1.0-self.u2) * self.r1 * self.n1) * self.n0
        return top / self.AvgFitness()
    
    def GetT21(self):
        top = (self.u1 * self.r0 * self.n0 + (1.0-self.u2) * self.r1 * self.n1) * self.n2
        return top / self.AvgFitness()
    
    def GetT02(self):
        top = (self.u2 * self.r1 * self.n1 + self.r2 * self.n2) * self.n0
        return top / self.AvgFitness()
        
    def GetT12(self):
        top = (self.u2 * self.r1 * self.n1 + self.r2 * self.n2) * self.n1
        return top / self.AvgFitness()
    
    def EventT10(self):
        self.n1 -= 1
        self.n0 += 1
        
    def EventT01(self):
        self.n1 += 1
        self.n0 -= 1
        
    def EventT20(self):
        self.n2 -= 1
        self.n0 += 1
        
    def EventT02(self):
        self.n2 += 1
        self.n0 -= 1
        
    def EventT12(self):
        self.n1 -= 1
        self.n2 += 1
        
    def EventT21(self):
        self.n1 += 1
        self.n2 -= 1
