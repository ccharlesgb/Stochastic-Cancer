# -*- coding: utf-8 -*-
"""
Created on Mon Oct 13 21:14:55 2014

@author: Jonny
"""

import scipy.integrate as integrate
import math
import matplotlib.pyplot as plt
import numpy as np
import Vi

datapointcount=100
N=10
r0=1.0
r1=1.01
r2=1.0
u1=0.1
u2=0.1
curpoint=100

rho1 = 0.0
rho3 = 0.0
V1 = 0.9
a = 0.0
b = 0.0
c = 0.0

def CalculateValues():
    global rho1
    global rho3
    global V1
    global a
    global b
    global c
    
    rho1 = 1.0 - ((r0 * (1.0-u1)) / (r1 + r0 * u1)) 
    rho1 = rho1 / (1.0 - math.pow((r0 * (1.0-u1) / (r1+r0 * u1)) , N))
    
    rho3 = 1.0 - (r0 * (1.0-u1) / (r2 + r1 * u2))
    rho3 = rho3 / (1.0 - math.pow((r0 * (1.0-u1) / (r2 + r1 * u2)) , N))
    
    
    
    V1=Vi.GetV_i(1, r0, r1, r2, u1, u2, N)
    
    a=N*u1*rho1
    b=N*u1*(1.0-V1-rho1)
    c=N*u2*(1.0-r1/r2)/(1.0-(math.pow((r1/r2) , N)))

CalculateValues()

print("a: {0}".format(a))
print("b: {0}".format(b))
print("c: {0}".format(c))
print("V1: {0}".format(V1))
print("rho1: {0}".format(rho1))

def gee(t):
    global c
    
    result = 1.0-math.exp(-(c*t))    
    
    return result    
    
    
    
def ell(t):
    global a
    global b
    global c
    
    result, error = integrate.quad(lambda z: gee(z)*a*math.exp(-(a+b)*(t-z)), 0, t)
    
    return result

def traj(t):
    result = (b/(a+b))*(1.0-math.exp(-(a+b)*t)) + ell(t)
    
    return result

minr1 = 0.2
maxr1 = 3.0
dataPointsY=[]
dataPointsT=[]
for i in range(0, datapointcount):
    
    r1 = ((float(i) / (datapointcount - 1)) * (maxr1 - minr1)) + minr1
    CalculateValues()
    dataPointsT.append(r1)
    dataPointsY.append(traj(100))

#
    
#t = np.arange(0.,  CurPoint., 1)
#y=traj(t)
plt.plot(dataPointsT, dataPointsY, linewidth=2.0)
plt.show
    