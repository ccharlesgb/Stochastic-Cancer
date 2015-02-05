# -*- coding: utf-8 -*-

import math
import matplotlib.pyplot as plt
import HaenoModel
'''
import SimTools

mySim = SimTools.Gillespie(10)
mySim.timeLimit = 1.0
mySim.r0 = 1.0
mySim.r1 = 1.0
mySim.r2 = 1.0

mySim.u1 = 0.1
mySim.u2 = 0.9
'''

myHaeno = HaenoModel.HaenoModel()
myHaeno.timeLimit = mySim.timeLimit

myHaeno.N = 10

myHaeno.r0 = 1.0
myHaeno.r1 = 0.5
myHaeno.r2 = 2.0

myHaeno.u1 = 0.1
myHaeno.u2 = 0.1

minr1 = 0.1
maxr1 = 3.0

dataPointsT = []
dataPointsB = []
dataPointsV = []
dataPointsR1 = []

datapointcount = 25
sdp = 1000

minr1 = 0.5
maxr1 = 3.0

for i in range(0, datapointcount):
    goodVi = 0
    
    myHaeno.r1 = ((float(i) / (datapointcount - 1)) * (maxr1 - minr1)) + minr1

    myHaeno.rho1 = myHaeno.Getrho1()
    myHaeno.rho3 = myHaeno.Getrho3()
    myHaeno.UpdateV() 
    
    b = myHaeno.N * myHaeno.u1 * (1.0 - myHaeno.V[1] - myHaeno.rho1)   
     
    dataPointsT.append(myHaeno.r1)
    dataPointsB.append(b)
    dataPointsR1.append(myHaeno.rho1)
        
    dataPointsV.append(myHaeno.V[1])

plt.subplot(311)
plt.plot(dataPointsT, dataPointsB, linewidth=4.0, label="b")
plt.title("B")
plt.subplot(312)    
plt.plot(dataPointsT, dataPointsV, linewidth=4.0, label="V1")   
plt.title("V1") 
plt.subplot(313)    
plt.plot(dataPointsT, dataPointsR1, linewidth=4.0, label="Rho1")    
plt.title("Rho1") 
        
        
        
        
#for si in range(0, sdp):
'''mySim.positive = 0
mySim.newLineage = 0
mySim.in0 = 10 - i
mySim.in1 = i
mySim.in2 = 0
mySim.N = 10

mySim.Simulate()

if mySim.positive == 1:
    goodVi += 1
elif mySim.newLineage == 0:
    goodVi += 1
    '''