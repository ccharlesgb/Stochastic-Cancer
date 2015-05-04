# -*- coding: utf-8 -*-
"""
Created on Mon May 04 12:49:13 2015

@author: Connor
"""

import math
import matplotlib.pyplot as plt
import SimUtil

dz = 0.001

minZ = 0.0
maxZ = 1.0

C = 1e-12

def GetD2Phi(phi):
    val = -math.sin(2.0 * phi)
    return val
    
z_vals = []
phi_vals = []
phi_prime_vals = []
phi_prime2_vals = []
phi_vals_new = []

phi_BEFORE = []

VAL_COUNT = int(maxZ / dz)

for i in range(0, VAL_COUNT):
    z_vals.append(SimUtil.SweepParameter(i,maxZ/dz,minZ,maxZ))
    phi_val = SimUtil.SweepParameterLog(i, maxZ/dz, 90.0, math.degrees(1.169))
    phi_vals.append(phi_val)
    phi_vals_new.append(phi_val)
    phi_BEFORE.append(phi_val)

plt.subplot(211)
plt.plot(z_vals, phi_vals)    

MAX_ITER = 1000

for i in range(0,MAX_ITER):
    #Update the d2phi

    #Update the phi
    for j in range(1, VAL_COUNT-1):
        cur_phi = phi_vals[j]
        fac = (dz * dz)/C
        phi_vals_new[j] = 0.5 * (fac * math.sin(2.0 * cur_phi) + phi_vals[j+1] + phi_vals[j-1])
        
    for j in range(1, VAL_COUNT-1):
        phi_vals[j] = phi_vals_new[j]

plt.subplot(212)
plt.plot(z_vals, phi_vals)

for i in range(0,VAL_COUNT):
    print(phi_BEFORE[i] - phi_vals[i])
        
    
