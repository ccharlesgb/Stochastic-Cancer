# -*- coding: utf-8 -*-
"""
Created on Mon Oct 13 21:14:55 2014

@author: Jonny
"""

import scipy.integrate as integrate
import math
import matplotlib.pyplot as plt
#import numpy as np
import Vi
import qk

datapointcount=15
N=10
r0=1.0
r1=0.3
r2=1.0
u1=0.1
u2=0.01
curpoint=100

t = 100.0

rho1 = 0.0
rho3 = 0.0
V1 = 0.0
a = 0.0
b = 0.0
c = 0.0

q0_points = 100
q0_vals = []

def CacheQ0():
    global q0_vals
    q0_vals = []
    for i in range(0, q0_points):
        curT = float(i)/(q0_points-1) * t
        q0_vals.append(qk.GetQ0(curT, r0, r1, r2, u1, u2, N))
        
def FetchQ0(time):
    index = (time / t) * (q0_points -1)
    #print("INDEX {0}".format(index))
    
    left_dex = int(math.floor(index))
    right_dex = int(math.ceil(index))
    left_val = q0_vals[left_dex]
    right_val = q0_vals[right_dex]
    
    frac = index - left_dex
    
    val = left_val + (right_val - left_val) * frac
    
    return val

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
    #print("V1: {0}".format(V1))
    
    a=N*u1*rho1
    b=N*u1*(1.0-V1-rho1)
    c=N*u2*(1.0-r1/r2)/(1.0-(math.pow((r1/r2) , N)))
    
    CacheQ0()

#CalculateValues()


def gee(t):
    global c
    #result = 1.0-math.exp(-(c*t))    
    #result = qk.GetQk(t, r0, r1, r2, u1, u2, N)
    result = FetchQ0(t)
    return result    
    
def ell(t):
    global a
    global b
    global c
    
    result, error = integrate.quad(lambda z: gee(z)*a*math.exp(-(a+b)*(t-z)), 0, t)
    return result

def traj(t):
    #print("A = {0} B = {1} B/A+B = {2}".format(a,b, b/(a+b)))
    result = (b/(a+b))*(1.0-math.exp(-(a+b)*t)) + ell(t)
    
    return result

minr1 = 0.1
maxr1 = 3.0
dataPointsY=[]
dataPointsT=[]
dataPoints1=[]
dataPoints2=[]

for i in range(0, datapointcount):
    r1 = ((float(i) / (datapointcount - 1)) * (maxr1 - minr1)) + minr1
    CalculateValues()
    dataPointsT.append(r1)
    dataPointsY.append(traj(t))
    #dataPoints1.append(traj(t) - ell(t))
    #dataPoints2.append(ell(t))
    #print(dataPointsY[i])
#

print("a: {0}".format(a))
print("b: {0}".format(b))
print("c: {0}".format(c))
print("V1: {0}".format(V1))
print("rho1: {0}".format(rho1))
#t = np.arange(0.,  CurPoint., 1)
#y=traj(t)
plt.plot(dataPointsT, dataPointsY, linewidth=4.0, label="X2(t)")
#plt.plot(dataPointsT, dataPoints1, linewidth=1.0, label = "1-exp()")
#plt.plot(dataPointsT, dataPoints2, linewidth=2.0, label = "L(t)")
plt.legend()
plt.show()

'''
ind=[] 
V=[]           
for i in range(0,N):
    ind.append(i)
    V.append(Vi.GetV_i(i, r0, r1, r2, u1, u2, N))
plt.plot(ind, V, linewidth=4.0, label="V_i")
'''