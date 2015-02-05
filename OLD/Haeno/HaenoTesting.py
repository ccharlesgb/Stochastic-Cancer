# -*- coding: utf-8 -*-

#Haeno Analysis Testing
#Sweeps r1

import math
import matplotlib.pyplot as plt
import HaenoModel
#import ./SimTools

myHaeno = HaenoModel.HaenoModel()

myHaeno.N = 10
myHaeno.r0 = 1.0
myHaeno.r1 = 1.0
myHaeno.r2 = 1.0

myHaeno.u1 = 0.1
myHaeno.u2 = 0.01

myHaeno.clampB = 1

datapointcount = 20

minr1 = 0.2
maxr1 = 3.0
dataPointsR1 = []
dataPointsX2 = []
dataPointsX0 = []
dataPointsX1 = []

for i in range(0, datapointcount):
    r1 = ((float(i) / (datapointcount - 1)) * (maxr1 - minr1)) + minr1
    dataPointsR1.append(r1)
    myHaeno.r1 = r1
    print("R0 = {0} R1 = {1}".format(myHaeno.r0,myHaeno.r1))
    X2 = myHaeno.GetX2()
    print("X2 is {0}".format(X2))
    dataPointsX2.append(X2)
    
    X0 = myHaeno.GetX0()
    dataPointsX0.append(X0)
    
    dataPointsX1.append(1.0 - X0 - X2)

#plt.plot(dataPointsR1, dataPointsX0, linewidth=1.0, label="X0(t)", linestyle = '--')
#plt.plot(dataPointsR1, dataPointsX1, linewidth=1.0, label="X1(t)", linestyle = '--')
plt.plot(dataPointsR1, dataPointsX2, linewidth=1.0, label="X2(t)", linestyle = '--')


