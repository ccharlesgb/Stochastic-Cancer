# -*- coding: utf-8 -*-

import scipy.integrate as integrate
import math
import matplotlib.pyplot as plt

class HaenoModel:
    def __init__(self):
        self.r0 = 1.0
        self.r1 = 1.1
        self.r2 = 1.0
        
        self.u1 = 0.1
        self.u2 = 0.1
        
        self.t = 100.0
        
        self.N = 100
        
        self.a = 0.0
        self.b = 0.0
        
        self.V = []
        
        self.Q0_points = 20
        self.Q0_vals = []
        
        self.rho1 = 0.0
        self.rh03 = 0.0
    
    def Geta(self):
        return self.N * self.u1 * self.rho1

    def Getb(self):
        return self.N * self.u1 * (1.0 - self.V[1] - self.rho1)

    def Getrho1(self):
        rho1 = 1.0 - ((self.r0 * (1.0-self.u1)) / (self.r1 + self.r0 * self.u1)) 
        rho1 = rho1 / (1.0 - math.pow((self.r0 * (1.0-self.u1) / (self.r1+self.r0 * self.u1)) , self.N))
        return rho1
        
    def Getrho3(self):
        rho3 = 1.0 - (self.r0 * (1.0-self.u1) / (self.r2 + self.r1 * self.u2))
        rho3 = rho3 / (1.0 - math.pow((self.r0 * (1.0-self.u1) / (self.r2 + self.r1 * self.u2)) , self.N))
        return rho3
        
    def UpdateV(self):
        self.V = []
        newV = []
        #Set up initial linear distribution of V
        for i in range(0,self.N+1):
            val = 1.0 - float(i)/self.N
            self.V.append(val)
            newV.append(val)
        
        MAX_ITER = 1000
        for itere in range(0,MAX_ITER):
            maxChange = -1.0
            for i in range(1, self.N):
                prv = self.V[i-1]
                nxt = self.V[i+1]
                Pi = float(i*(self.N-i))
                Pi /= (i*self.r1*(1.0-self.u2) + self.r0*(self.N-i))
                
                new = self.r1 * (1.0 - self.u2) * nxt + self.r0 * prv
                new = new / (self.r1 * self.u2 * self.rho3/Pi + self.r1 * (1.0 - self.u2) + self.r0)
                
                if abs(new - self.V[i]) > maxChange:
                    maxChange = abs(new-self.V[i])                
                
                newV[i] = new
                #self.V[i]= new
            
            for i in range(1, self.N):          
                self.V[i] = newV[i]
                newV[i] = self.V[i]
                
            if maxChange < 1e-15:
                break
        
    def CalculateParameters(self):
        self.rho1 = self.Getrho1()
        self.rho3 = self.Getrho3()
        
        self.UpdateV()
        self.UpdateQ0()
        
        self.a = self.Geta()
        self.b = self.Getb()
        if self.b < 0.0:
            self.b = 0.0
        
        print("A = {0} B = {1}".format(self.a, self.b))
    
    def GetLambdaK(self, k):
        res = (self.N - k) * ((self.N - k) * self.r1 * self.u2 + self.r2 * k)
        res /= (self.N-k)*self.r1 + self.r2 * k
        return res
    
    
    def GetMuK(self, k):
        res = k * (self.N - k) * self.r1 * (1.0 - self.u2)
        res /= (self.N - k) * self.r1 + self.r2 * k
        return res    
    
    def UpdateQ0(self):
        self.Q0_vals = []
        reachedMax = 0 #If a q0 = 1.0 it will remain 1.0
        for t_it in range(0, self.Q0_points):

            if reachedMax == 1:
                self.Q0_vals.append(1.0)
                continue
            
            maxT = float(t_it) / (self.Q0_points - 1) * self.t
    
            #Initialise Q with I conditions
            q = []
            q_dots = []
            for i in range(0,self.N+1):
                q.append(0.0)
                q_dots.append(0.0)
            q[self.N] = 1.0
    
            curT = 0.0
            deltaT = 0.01    
            while curT < maxT:
                for k in range(0,self.N):
                    lamk = self.GetLambdaK(k)
                    muk = self.GetMuK(k)
                    if k > 0:
                        q_dots[k] = lamk * q[k+1] - (lamk + muk) * q[k] + muk * q[k-1]
                    else:
                        q_dots[k] = lamk * q[k+1] - (lamk + muk) * q[k]
                for k in range(0,self.N):
                    newQ = q[k] + q_dots[k] * deltaT
                    q[k] = newQ
                curT += deltaT
            
            if (abs(q[0] - 1.0) < 0.0000001):
                reachedMax = 1 #Stop calculating all other q0 will be 1.0 for extra time
                print("REACHED MAX")
            
            print("Q_{0} = {1}".format(maxT, q[0]))
            self.Q0_vals.append(q[0])
    
    #Fetch Cached Q0
    def Q0(self, t):
        index = (t / self.t) * (self.Q0_points - 1)
        #print("INDEX {0}".format(index))
        
        left_dex = int(math.floor(index))
        right_dex = int(math.ceil(index))
        left_val = self.Q0_vals[left_dex]
        right_val = self.Q0_vals[right_dex]
        
        frac = index - left_dex
        
        #Linerly interlopate between Q0 vals
        val = left_val + (right_val - left_val) * frac
        return val

    def GetTunnelTerm(self):
        res = (self.b / (self.a + self.b)) * (1.0 - math.exp(-(self.a + self.b) * self.t))
        return res
        
    def GetL(self):
        result, error = integrate.quad(lambda z: self.Q0(z) * self.a * math.exp(-(self.a + self.b)*(self.t - z)), 0, self.t)
        return result
    
    def GetX2(self):
        self.CalculateParameters()
        tunnelTerm = self.GetTunnelTerm()
        L_term = self.GetL()
        print("TT {0} LT {1}".format(tunnelTerm, L_term))
        return tunnelTerm + L_term
        

