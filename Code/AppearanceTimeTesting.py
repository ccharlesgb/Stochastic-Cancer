# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 16:03:47 2015

@author: Connor
"""

import wright_fisher
import matplotlib.pyplot as plt
import math
import SimUtil
import TauSolver
import TauLeapParam
import MatTools
import random
    
cellTypes = 21
population = 1e9

myWF = wright_fisher.wright_fisher(cellTypes)
myHist = TauLeapParam.Hist(cellTypes)
myParam = TauLeapParam.Params(cellTypes)

myParam.d = 100

myWF.stopAtAppear = 1

myParam.n0[0] = population
myWF.history = myHist
myWF.timeLimit = 1000000
myWF.useApproxTheta = 0
myWF.params = myParam

SPD = 1
DPC = 1

myHist.SPD = SPD

s = 1e-2
myParam.SetCompoundFitness(s)   
myParam.SetUAll(1e-4)

#Sinusoidal mutation rate
f = 0.25
for i in range(1, cellTypes):
    u_log = -7.0 + 4.0*math.cos(i*f * 2.0 * math.pi)
    u = math.pow(10.0, u_log)
    print(u_log),
    myParam.SetU(i,u)
print("")

'''
#Exponential decrease mutation rate
for i in range(1, cellTypes):
    u =  SimUtil.SweepParameterLog(i-1,cellTypes-1,1e-6,1e-9)
    u_log = math.log10(u)
    print(u_log),
    myParam.SetU(i,u)
'''
'''
#Exponential increase mutation rate
for i in range(1, cellTypes):
    u =  SimUtil.SweepParameterLog(i-1,cellTypes-1,1e-9,1e-6)
    u_log = math.log10(u)
    print(u_log),
    myParam.SetU(i,u)
'''

#Switch
u_before = 1e-6
u_after = 1e-9
switch = 10
for i in range(1, cellTypes):
    if (switch + 1) > i:
        myParam.SetU(i, u_before)
    else:
        myParam.SetU(i, u_after)


mySolver = TauSolver.Solver(myParam)
mySolver.CacheX0()

avgAppear = [0.0] * cellTypes
theoreticalAppear = [0.0] * cellTypes
theoreticalAppear_NEW = [0.0] * cellTypes
theoreticalAppear2 = [0.0] * cellTypes
theoreticalAppear_MODEL = [0.0] * cellTypes

Eq12 = []
Eq12_Old = []
for i in range(0, cellTypes):
    Eq12.append([])
    Eq12_Old.append([])

mySolver.tau_hist = [0.0]*(cellTypes-1)
for curPoint in range(0,SPD):
    myWF.Simulate()
    print("Sim Done")
    last_t = 0.0
    #Loop through each step
    for i in range(0, cellTypes):
        found = False

        for t in range(len(myHist.histArray[i])):
            if found:
                continue
            count = myHist.histArray[i][t]
            if count > 0:
                if curPoint == 0 and i < cellTypes -1:
                    mySolver.tau_hist[i] = t - last_t
                    last_t = t
                    for t_prime in range(t, t + 300):
                        
                        val = (mySolver.GetXJ(t_prime-t,i) * myParam.N)
                        
                        Eq12[i].append(val)
                        
                        gam = math.sqrt(2.0 * math.log(myParam.N))
                        Eq12_Old[i].append(mySolver.GetXJ_NEW(i, t_prime-t) * myParam.N)
                        #Clamp to N
                        if Eq12_Old[i][len(Eq12_Old[i]) - 1 ] > myParam.N:
                            Eq12_Old[i][len(Eq12_Old[i])-1] = 0.0
                        if Eq12[i][len(Eq12[i]) -1 ] > myParam.N:
                            Eq12[i][len(Eq12[i])-1] = 0.0
                avgAppear[i] += t
                found = True

theory_NEW_total = 0.0
theory_model_total = 0.0
#VERY IMPORTANT THESE ARE CHECKED
mySolver.tau_hist = [0.0]*(cellTypes-1)
for i in range(0, cellTypes - 1): #From 1 to 20
    theoreticalAppear[i+1] = (i+1) * mySolver.GetTauOriginal()
    theoreticalAppear_NEW[i+1] = theory_NEW_total + mySolver.GetTauCorrect(i)
    theory_NEW_total = theoreticalAppear_NEW[i+1]
    theoreticalAppear_MODEL[i+1] = theory_model_total + mySolver.GetTauRecursive(i)
    theory_model_total = theoreticalAppear_MODEL[i+1]
    
for i in range(0, cellTypes):
    avgAppear[i] /= SPD

print("ORIG", mySolver.GetWaitingTimeOriginal(cellTypes-1))

j_i = int(round(-math.log(myParam.N)/math.log(myParam.GetU(1)*myParam.d)))
frac_part = j_i - math.floor(j_i)
first_whole = int(math.ceil(j_i))

#Take the fractioanal part of the one we are missing
frac_tau = frac_part * mySolver.GetTauNeglect(int(math.floor(j_i)))

total = frac_tau
theoreticalAppear2[first_whole] = frac_tau
for i in range(first_whole, cellTypes-1):
    total = total + mySolver.GetTauNeglect(i)
    theoreticalAppear2[i+1] = total

y = range(0,cellTypes)
plt.figure()
plt.subplot(211)
plt.plot(avgAppear, 'o', label = "Simulation")
#plt.plot(theoreticalAppear, label = "Original Method")
plt.plot(theoreticalAppear2, '--', label = "Neglecting", linewidth = 2)
plt.plot(theoreticalAppear_NEW,label="Single Correction", linewidth = 2)
plt.plot(theoreticalAppear_MODEL, ':', label="Recursive", linewidth = 2)
plt.xlabel("i")
plt.ylabel("Appearance Time")
plt.legend(loc = 2)
plt.title("U = {0}, {1}, {2},{3}...".format(myParam.u[0],myParam.u[1],myParam.u[2],myParam.u[3]))

plt.subplot(212)
for i in range(0, cellTypes-1):
   plt.plot(myHist.tTotal, myHist.histArrayTotal[i])
   #plt.plot(range(int(avgAppear[i]),int(avgAppear[i]) + 300), Eq12[i], ':', linewidth = 2.0)
   #plt.plot(range(int(avgAppear[i]),int(avgAppear[i]) + 300), Eq12_Old[i], '--')
#plt.yscale("log")
plt.show()

data = dict()
data["appear_sim"] = avgAppear
data["appear_orig"] = theoreticalAppear
data["appear_correct"] = theoreticalAppear_NEW
data["appear_neglect"] = theoreticalAppear2
data["appear_recursive"] = theoreticalAppear_MODEL

MatTools.SaveDict2(data, prefix = "VARUSWITCH", spd = SPD, dpc = DPC, params = myParam.GetFileString())



