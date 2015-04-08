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

SPD = 5
DPC = 1

s = 1e-1
for i in range(0,cellTypes):
    myParam.r[i] = math.pow(1.0 + s, i)
    
myParam.u = [1e-8] * cellTypes

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

for curPoint in range(0,SPD):
    myWF.Simulate()
    print("Sim Done")
   
    #Loop through each step
    for i in range(0, cellTypes):
        found = False
        for t in range(len(myHist.histArray[i])):
            if found:
                continue
            count = myHist.histArray[i][t]
            if count > 0:
                if curPoint == 0:
                    for t_prime in range(t, t + 300):
                        
                        val = mySolver.GetXJ(t_prime-t,i) * myParam.N    
                        
                        Eq12[i].append(val)
                        
                        gam = math.sqrt(2.0 * math.log(myParam.N))
                        Eq12_Old[i].append(math.exp((t_prime - t) * s * gam))
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
for i in range(0, cellTypes - 1): #From 1 to 20
    theoreticalAppear[i+1] = (i+1) * mySolver.GetTauOriginal()
    theoreticalAppear_NEW[i+1] = theory_NEW_total + mySolver.GetTau(i)
    theory_NEW_total = theoreticalAppear_NEW[i+1]
    
    theoreticalAppear_MODEL[i+1] = theory_model_total + mySolver.GetTauModel(i)
    theory_model_total = theoreticalAppear_MODEL[i+1]
    
for i in range(0, cellTypes):
    avgAppear[i] /= SPD

print("ORIG", mySolver.GetWaitingTimeOriginal(cellTypes-1))

j_i = int(round(-math.log(myParam.N)/math.log(myParam.u[0]*myParam.d)))
for i in range(j_i, cellTypes):
    theoreticalAppear2[i] = (i-j_i) * mySolver.GetTauNeglect()

plt.figure()
plt.subplot(211)
plt.plot(avgAppear, 'o')
plt.plot(theoreticalAppear)
plt.plot(theoreticalAppear2, '--')
plt.plot(theoreticalAppear_NEW, '^')
plt.plot(theoreticalAppear_MODEL, '*')
plt.xlabel("i")
plt.ylabel("Appearance Time")

plt.subplot(212)
for i in range(0, cellTypes):
   plt.plot(myHist.tHist[0::1], myHist.histArray[i][0::1])
   plt.plot(range(int(avgAppear[i]),int(avgAppear[i]) + 300), Eq12[i], ':', linewidth = 2.0)
   plt.plot(range(int(avgAppear[i]),int(avgAppear[i]) + 300), Eq12_Old[i], '--')
plt.yscale("log")

plt.yscale("log")
plt.show()

data = dict()
data["appear_sim"] = avgAppear
data["appear_orig"] = theoreticalAppear
data["appear_transient"] = theoreticalAppear_NEW
data["appear_neglect"] = theoreticalAppear2

MatTools.SaveDict2(data, spd = SPD, dpc = DPC, params = myParam.GetFileString())



