# -*- coding: utf-8 -*-
"""
Created on Mon Nov 24 22:14:28 2014
@author: Jonny
"""
import numpy as np
import math

class systematic:
    def __init__(self, sim):
        self.sim=sim
        self.reset()
        self.threshold = 1e-6
        self.deltaT = 0.1        
        
    def reset(self): #function to reset system parameters after a run
         self.w = np.zeros((self.sim.N +1, self.sim.N + 1 ))
         self.w_dots = np.zeros((self.sim.N +1, self.sim.N + 1 ))
         self.w_dots_cached = []
         
         self.curT = 0.0
         self.tmax = self.sim.timeLimit
         self.tConverge = self.tmax
         
         self.integral = 0.0
    
    def cumulative(self): #function to 
        self.w[0][self.sim.N] = 1.0

        reachedMax = 0
        while self.curT < self.tmax: #loop to save time if the system has converged
            if self.sim.preSim != 0:
                self.sim.preSim(self)
            
            if reachedMax == 1:
              for i in range(0,self.sim.N+1):
                 for j in range(0,self.sim.N+1):
                     self.w[i][j] = 1.0
                     self.w_dots[i][j] = 0.0
              self.w_dots_cached.append(0.0)
              self.curT += self.deltaT
              continue
            
            self.Update_W_dot()
            for i in range(0,self.sim.N+1):
               for j in range(0,self.sim.N+1):
                   self.w[i][j] = self.w[i][j] + self.deltaT * self.w_dots[i][j]
           
            w_dot_term=self.w_dots[0][0]
            print(w_dot_term)
            self.w_dots_cached.append(w_dot_term)
            
            self.integral += w_dot_term * self.curT * self.deltaT           
            
            self.curT += self.deltaT
            
            if (abs(self.w[0][0] - 1.0) < self.threshold) and reachedMax == 0:
                reachedMax = 1
                print("System has converged!") 
                self.tConverge = self.curT
                #self.integral += (self.tmax - self.curT) * 
                break
                                

    def Get_fix_time(self):
        self.reset()
        print("Cumulative running..")
        self.cumulative()
        return self.integral

    def Get_w_dot_00(self,t):
        if (t > self.tConverge):
            return 0.0
        index = (t / self.tConverge) * (len(self.w_dots_cached) - 1.0)
        #print("INDEX {0}".format(index))
        
        left_dex = int(math.floor(index))
        right_dex = int(math.ceil(index))
        left_val = self.w_dots_cached[left_dex]
        right_val = self.w_dots_cached[right_dex]
        
        frac = index - left_dex
        
        #Linerly interlopate between Q0 vals
        val = left_val + (right_val - left_val) * frac
        return val
        
    def Update_W_dot(self):
        for i in range(0,self.sim.N+1):
            #self.i=i
            for j in range(0,self.sim.N+1):
                #self.j=j                
                if i+j>self.sim.N:
                    continue
                self.sim.n0 = self.sim.N - i - j
                self.sim.n1 = i
                self.sim.n2 = j
                
                self.w_dots[i][j] = self.sim.GetT01()  * (self.Get_w_ij(i+1, j) - self.Get_w_ij(i,j) )
                self.w_dots[i][j]+= self.sim.GetT10() * (self.Get_w_ij(i-1, j) - self.Get_w_ij(i,j) )
                self.w_dots[i][j]+= self.sim.GetT02()  * (self.Get_w_ij(i  ,j+1) - self.Get_w_ij(i,j) )
                self.w_dots[i][j]+= self.sim.GetT20() * (self.Get_w_ij(i  ,j-1) - self.Get_w_ij(i,j) )
                self.w_dots[i][j]+= self.sim.GetT21() * (self.Get_w_ij(i+1,j-1) - self.Get_w_ij(i,j) )
                self.w_dots[i][j]+= self.sim.GetT12() * (self.Get_w_ij(i-1,j+1) - self.Get_w_ij(i,j) )
        
    def Get_w_ij(self,i,j):
        if i+j > self.sim.N:
            return 0.0
        if i < 0:
            return 0.0
        if j < 0:
            return 0.0
        element = self.w[i][j]
        return element