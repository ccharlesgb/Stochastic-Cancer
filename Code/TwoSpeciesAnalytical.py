# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 11:58:48 2015

@author: Jonny
"""
import numpy as np

class TwoSpeciesAnalytical:
    def __init__(self, params):
        self.sim = 0
        self.params = params
        
        self.cached_gamma = []
        self.cached_inverse_TJ_Plus = []

        self.N = 100         
        self.r0 = 1.0
        self.r1 = 1.0
        self.u1 = 0.1        
        
        self.j = 0

    def reset(self):
        self.N = self.params.n0[0]
        self.r0 = self.params.r[0]
        self.r1 = self.params.r[1]
        self.u1 = self.params.GetU(1)
        
        self.cached_gamma = [0]*(self.N)
        self.cached_inverse_TJ_Plus = [0]*(self.N)
        self.cached_gamma_product = np.zeros((self.N + 1, self.N +1))
        
          
    def cache_gamma(self):
        for j in range(0, self.N):            
            self.cached_gamma[j] = self.GetGammaJ(j)
        
    def cache_inverse_TJ_Plus(self):
        for j in range(0,self.N):
            self.cached_inverse_TJ_Plus[j] = 1.0/self.GetTJplus(j)
     
    def cache_gamma_prod(self):
        #first ininitalise all terms where l is less than k to be 1 - the empty product        
        for m in range(0,self.N):   
            print("Caching the gamma product.. point {0}/{1}".format(m,self.N))            
            for k in range(0, m):
                self.cached_gamma_product[m][k] = 1.0
            #now initialise the first valid product
            self.cached_gamma_product[m][m] = self.cached_gamma[m]
            #now loop through endpoints
            for k in range(m+1, self.N):
                self.cached_gamma_product[m][k] = self.cached_gamma_product[m][k-1]*self.cached_gamma[k]
        for k in range(0,self.N):
            self.cached_gamma_product[self.N][k] = 1.0

    def GetCachedGammaProd(self, l, k):
        if(l>k):
            return 1.0
        if(self.cached_gamma_product[l-1] == 0.0):
            print("Warning - camma product is very close to zero")
            return 0.0
        return self.cached_gamma_product[k]/self.cached_gamma_product[l-1]
     
    def GetAverageFitness(self, j):
        return ( (self.N - j)*self.r0 + j*self.r1)
        
    def GetTJplus(self, j):
        top = (self.N-j)*( (j * self.r1) + self.u1 * self.r0 * (self.N - j))
        #return top / (self.params.GetAvgFit([self.N, self.N-j]) )
        return top / self.GetAverageFitness(j)
    
    #probability for j -> j-1
    def GetTJminus(self, j):
        top = j * (self.N - j) * self.r0 * (1.0 - self.u1)
        #return top / (self.params.GetAvgFit([self.N, self.N-j]) )
        return top / self.GetAverageFitness(j)
    
    
    def GetGammaJ(self, j):
        return self.GetTJminus(j)/ self.GetTJplus(j)
        
    def getFixTime(self,i):
        print("u1 is: {0}, r0 is: {1}, r1 is: {2}".format(self.u1, self.r0, self.r1))        
        
        self.reset()        
        print("Caching gamma...")        
        self.cache_gamma()      
        #print(self.cached_gamma)
        print("Caching TJ Plus...")
        self.cache_inverse_TJ_Plus()
        print("Caching gamma product")
        self.cache_gamma_prod()
        print(self.cached_gamma_product)        

        sum_k = 0     
        for k in range(i,self.N):
           print("k is: {0}".format(k))
           first_term = (1.0/(self.N * self.u1) ) * self.cached_gamma_product[1][k] 
           
           second_term = 0         
           for m in range(1, k+1):        
               second_term += self.cached_inverse_TJ_Plus[m]*self.cached_gamma_product[m+1][k]
                
           sum_k += first_term + second_term
        
        return sum_k
'''
def getFixTimeNew(sim,i):

    sum_k = 0     
    for k in range(i,sim.N):
        #print("k is: {0}".format(k))
        first_term = 1.0/(sim.N * sim.u1)    
        second_term = 0         
        for m in range(1, k+1):        
            sim.j = m            
            first_term *= GetGammaj(sim)
        
            internal_product = 1.0/sim.GetTJplus()
            for l in range(m+1, k+1):
                sim.j = l
                internal_product*= GetGammaj(sim)
            second_term += internal_product
            
        sum_k += first_term + second_term
        
    return sum_k


'''
    