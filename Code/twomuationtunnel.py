# -*- coding: utf-8 -*-
"""
Created on Mon Oct 13 21:14:55 2014

@author: Jonny
"""

import scipy.integrate as integrate
import math
import matplotlib.pyplot as plt
import numpy as np

datapointcount=100
N=100
r0=1.0
r1=1.1
r2=1.0
u1=0.1
u2=0.01
V1=0.7
curpoint=100

rho1=(1.0-(r0*(1.0-u1)/(r1+r0*u1)))/(1.0-math.pow((r0*(1.0-u1)/(r1+r0*u1)) , N))
rho3=(1.0-(r0*(1.0-u1)/(r2+r1*u2)))/(1.0-math.pow((r0*(1.0-u1)/(r2+r1*u2)) , N))

V = []
V_lin = []

for i in range(0,N):
    V.append(1.0 - float(i+1) / N)
    V_lin.append(V[i])


MAX_ITER = 1000
MIN_CHANGE = 1e-6

V1_old = V[1]

for it in range(0, MAX_ITER):
    for i in range(1, N-1):
        nxt = V[i+1]
        prv = V[i-1]
        P_i = float(i*(N-i)) / (i*r1*(1.0-u2) + r0*(N-i))
        new = (r1*(1-u2) * nxt + r0 * prv)
        new = new / ((r1*u2*rho3)/P_i + r1*(1-u2) + r0)
        V[i] = new
    
    if abs(V[1] - V1_old) < MIN_CHANGE:
        print("Converged in {0} steps.".format(i))
        break
    V1_old = V[1]

V1=V[1]

a=N*u1*rho1
b=N*u1*(1.0-V1-rho1)
c=N*u2*(1.0-r1/r2)/(1.0-(math.pow((r1/r2) , N)))



print("a: {0}".format(a))
print("b: {0}".format(b))
print("c: {0}".format(c))
print("rho1: {0}".format(rho1))

def gee(t):
    global c
    
    result = 1.0-math.exp(-(c*t))    
    
    return result    
    
    
    
def ell(t):
    global a
    global b
    global c
    
    result = integrate.quad(lambda z: gee(z)*a*math.exp(-(a+b)*(t-z)), 0, t)
    
    return result[0]

def traj(t):
    result = (b/(a+b))*(1.0-math.exp(-(a+b)*t)) + ell(t)
    
    return result

dataPointsY=[]
dataPointsT=[]
for i in range(0, datapointcount):
    dataPointsT.insert(i,i)
    dataPointsY.insert(i,traj(i))


#
    
#t = np.arange(0.,  CurPoint., 1)
#y=traj(t)
plt.plot(dataPointsT, dataPointsY, linewidth=2.0)
plt.show
    