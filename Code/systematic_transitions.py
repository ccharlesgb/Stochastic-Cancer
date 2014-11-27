# -*- coding: utf-8 -*-
"""
Created on Mon Nov 24 22:14:28 2014

@author: Jonny
"""
import numpy as np
import scipy.integrate as integrate
import math

class systematic:
    def __init__(self, sim):
        self.reset(sim)
    
    def reset(self,sim):
         self.r0 = sim.r0
         self.r1 = sim.r1
         self.r2 = sim.r2
         self.u1 = sim.u1
         self.u2 = sim.u2
         
         self.N = sim.in0
         self.i = 0
         self.j = 0
         self.w = np.zeros((self.N +1, self.N + 1 ))
         self.w_dots= np.zeros((self.N +1, self.N + 1 ))
         self.w_dots_cached=[]
         
         self.tmax = sim.timeLimit
         self.threshold = 1e-6
    
    
    
    def cumulative(self):
        #self.Q0_vals = []
        #reachedMax = 0 #If a q0 = 1.0 it will remain 1.0
        
        #Initialise W with Initial conditions
        #self.w = np.zeros((self.N +1, self.N + 1 ))
        self.w[0][self.N] = 1.0
        
        #self.w_dots = np.zeros((self.N + 1, self.N + 1))
        
        curT = 0.0
        deltaT = 0.1
        reachedMax = 0
        while curT < self.tmax: #loop to save time if the system has converged
            if reachedMax == 1:
              for i in range(0,self.N+1):
                 for j in range(0,self.N+1):
                     self.w[i][j] = 1.0
                     self.w_dots[i][j] = 0.0
              self.w_dots_cached.append(0.0)
              curT += deltaT
              continue
            
            
            self.Update_W_dot()
            for i in range(0,self.N+1):
               for j in range(0,self.N+1):
                   self.w[i][j] = self.w[i][j] + deltaT * self.w_dots[i][j]
            curT += deltaT
           
            w_dot_term=self.w_dots[0][0]
            #print(self.w)
            self.w_dots_cached.append(w_dot_term)
            if (abs(self.w[0][0] - 1.0) < self.threshold) and reachedMax == 0:
                reachedMax = 1
                print("System has converged!")
                                

    def Get_fix_time(self,sim):
        self.reset(sim) 
        print("Cumulative running..")
        self.cumulative()
        print("Integrating...")
        result, error = integrate.quad(lambda t: self.Get_w_dot_00(t)*t, 0, self.tmax)
        
        deltaT = 0.1
        curT = 0.0

        res = 0.0
        
        while(curT < self.tmax):
            change = self.Get_w_dot_00(curT) * curT * deltaT
            res += change
            curT += deltaT
            
        print("RESULT: {0} OR {1}".format(result, res))
        
        return result

    def Get_w_dot_00(self,t):
        index = (t / self.tmax) * (len(self.w_dots_cached) - 1.0)
        #print("INDEX {0}".format(index))
        
        left_dex = int(math.floor(index))
        right_dex = int(math.ceil(index))
        left_val = self.w_dots_cached[left_dex]
        right_val = self.w_dots_cached[right_dex]
        
        frac = index - left_dex
        
        #Linerly interlopate between Q0 vals
        val = left_val + (right_val - left_val) * frac
        return val
        
    def getB(self):
        fitness = self.r0*(self.N-self.i-self.j) + self.r1*self.i +self.r2*self.j
        #print("i {0} {1} {2}".format(self.N-self.i-self.j, self.i, self.j))
        #print("R0 {0}".format(self.r0))
        #print("R2 {0}".format(self.r2))
        #print("R1 {0}".format(self.r1))
        return 1.0/fitness
        
    def Get_i_plus_j_stay(self):
        births = self.i*self.r1*(1.0-self.u2) + (self.N-self.i-self.j)*self.r0*self.u1       
        death = (self.N-self.i-self.j)*self.getB()
        return births*death
        
    def Get_i_minus_j_stay(self):
        births = self.r0*(1.0-self.u1)
        death = self.i*self.getB()
        return births*death
        
    def Get_i_stay_j_plus(self):
        births = self.j*self.r2 + self.i*self.r1*self.u2
        death = (self.N-self.i-self.j)*self.getB()
        return births*death
        
    def Get_i_stay_j_minus(self):
        births = (self.N-self.i-self.j)*self.r0*(1.0-self.u1)
        death = self.j*self.getB()
        return births*death
        
    def Get_i_plus_j_minus(self):
        births = self.i*self.r1*(1.0-self.u2) + (self.N-self.i-self.j)*self.r0*self.u1
        death = self.j*self.getB()
        return births*death
        
    def Get_i_minus_j_plus(self):        
        births = self.j*self.r2 + self.i*self.r1*self.u2
        death = self.i*self.getB()
        return births*death
        
    def Update_W_dot(self):
        for i in range(0,self.N+1):
            self.i=i
            for j in range(0,self.N+1):
                self.j=j                
                if i+j>self.N:
                    continue
                
                self.w_dots[i][j] = self.Get_i_plus_j_stay()  * (self.Get_w_ij(i+1, j) - self.Get_w_ij(i,j) )
                self.w_dots[i][j]+= self.Get_i_minus_j_stay() * (self.Get_w_ij(i-1, j) - self.Get_w_ij(i,j) )
                self.w_dots[i][j]+= self.Get_i_stay_j_plus()  * (self.Get_w_ij(i  ,j+1) - self.Get_w_ij(i,j) )
                self.w_dots[i][j]+= self.Get_i_stay_j_minus() * (self.Get_w_ij(i  ,j-1) - self.Get_w_ij(i,j) )
                self.w_dots[i][j]+= self.Get_i_plus_j_minus() * (self.Get_w_ij(i+1,j-1) - self.Get_w_ij(i,j) )
                self.w_dots[i][j]+= self.Get_i_minus_j_plus() * (self.Get_w_ij(i-1,j+1) - self.Get_w_ij(i,j) )
        
        
    def Get_w_ij(self,i,j):
        if i+j > self.N:
            return 0.0
        if i < 0:
            return 0.0
        if j < 0:
            return 0.0
        element = self.w[i][j]
        return element
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        