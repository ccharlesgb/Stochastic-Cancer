# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 18:41:05 2015

@author: Connor
"""

import math
import matplotlib.pyplot as plt

d = 100
s = 0.01
u = 1e-6
N = 1e9

print (u)

def F_Tau(tau):
    res = (N * u * d)/s
    res *= (float(u * d)/s * (math.exp(s*tau) - 1.0 ) - u * d * tau)
    return res
    
def F_Tau_orig(tau):
    res = (u * d) / s * math.exp(s * tau) - 1
    return res
    
dataX = []
dataY = []
dataY_orig = []

DC = 100
minTau = 0
maxTau = 10

onePoint = 0.0

for i in range(0, DC):
    tau = float(i) / (DC - 1) * (maxTau - minTau) + minTau
    dataX.append(tau)
    dataY.append(F_Tau(tau))
    dataY_orig.append(F_Tau_orig(tau))
    print("Tau = {0} Data = {1}".format(tau, dataY[i]))
    if abs(dataY[i] - 1.0) < 0.01 or dataY[i] > 1.0:
        onePoint = tau
        print("FOUND ONE POINT = {0}".format(onePoint))
        break
    
plt.figure()
plt.plot(dataX,dataY)
#plt.plot(dataX,dataY_orig)
    
    
    
    
    
    