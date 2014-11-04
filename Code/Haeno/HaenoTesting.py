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
myHaeno.r1 = 1.1
myHaeno.r2 = 1.0

myHaeno.u1 = 0.1
myHaeno.u2 = 0.1

datapointcount = 40

minr1 = 0.2
maxr1 = 3.0
dataPointsY=[]
dataPointsT=[]

for i in range(0, datapointcount):
    r1 = ((float(i) / (datapointcount - 1)) * (maxr1 - minr1)) + minr1
    dataPointsT.append(r1)
    myHaeno.r1 = r1
    print("R0 = {0} R1 = {1}".format(myHaeno.r0,myHaeno.r1))
    X2 = myHaeno.GetX2()
    print("X2 is {0}".format(X2))
    dataPointsY.append(X2)


plt.plot(dataPointsT, dataPointsY, linewidth=1.0, label="X2(t)", linestyle = '--')