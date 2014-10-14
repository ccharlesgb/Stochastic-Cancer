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
N=100
r0=1.0
r1=1.1
r2=1.0
u1=0.1
u2=0.01
curpoint=100

rho1=(1.0-(r0*(1.0-u1)/(r1+r0*u1)))/(1.0-math.pow((r0*(1.0-u1)/(r1+r0*u1)) , N))
rho3=(1.0-(r0*(1.0-u1)/(r2+r1*u2)))/(1.0-math.pow((r0*(1.0-u1)/(r2+r1*u2)) , N))



V1=Vi.GetV_i(1, r0, r1, r2, u1, u2, N)

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
    