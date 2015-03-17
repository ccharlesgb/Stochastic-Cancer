# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 16:03:47 2015

@author: Connor
"""

import wright_fisher
import matplotlib.pyplot as plt
import math
import SimUtil
import NumericalTau_TEST

def GetTau2(param):
    s = param.r[1] - param.r[0]
    logs = math.log(1.0 + s / (param.u[0] * param.d)*math.sqrt(2.0*math.log(param.popSize)))
    top = math.pow(logs,2.0)
    bottom = 2.0 * s * math.log(param.popSize)
    return top/bottom

def GetTau(param):
    s = param.r[1] - param.r[0]
    top = math.pow(math.log(s / (param.u[0] * param.d)),2.0)
    bottom = 2.0 * s * math.log(param.popSize)
    return top/bottom
    
cellTypes = 21
population = 1e9

myWF = wright_fisher.wright_fisher()
myHist = wright_fisher.wf_hist(cellTypes)
myParam = wright_fisher.wright_fisher_params(cellTypes)

myParam.d = 100

myWF.stopAtAppear = 1

myParam.iN[0] = population
myWF.history = myHist
myWF.stepLimit = 1000000
myWF.useApproxTheta = 0
myWF.params = myParam

myParam.u[0] = 1e-7

SPD = 1
DPC = 5

s = 1e-4
for i in range(0,cellTypes):
    myParam.r[i] = math.pow(1.0 + s, i)

myWF.reset()

avgAppear = [0.0] * cellTypes
theoreticalAppear = [0.0] * cellTypes
theoreticalAppear2 = [0.0] * cellTypes
theoreticalAppear_NEW = [0.0] * cellTypes

Eq12 = []
Eq12_Old = []
for i in range(0, cellTypes):
    Eq12.append([])
    Eq12_Old.append([])

for curPoint in range(0,SPD):
    myWF.reset()
    myWF.Simulate()       
   
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
                        x_j_1 = 0.0
                        if i > 0:
                            x_j_1 = float(myHist.histArray[i-1][t_prime-t]) / myParam.popSize
                            start = 1.0
                            if i == 1:
                                start = 1e4/myParam.popSize
                            u = myParam.u[0]
                            d = myParam.d
                            s = 0.01
                            T = t_prime - t
                            x_j = float(myHist.histArray[i][T])/myParam.popSize
                            jbar = myHist.avgJHist[t_prime - t]
                            if i == 1:
                                x_j_1 = 1.0
                            if i == 2:
                                x_j_1 = 1e4 / myParam.popSize
                            if i == 3:
                                x_j_1 = 1 / myParam.popSize
                            gam = math.sqrt(2.0 * math.log(myParam.popSize))    
                            #val = ((u * d * x_j_1 + (s * i * gam + u*d) * start/1e9) * math.exp((s-start*u*d) * gam * i * T) - u * d * x_j_1) / (s * i * gam + u*d)    
                            val = NumericalTau_TEST.GetXJ(T,i,myParam)
                            val = val * myParam.popSize                 
                        else:
                            val = count * math.exp((s + myParam.u[0] * myParam.d * x_j_1)*(t_prime-t) * (i - myHist.avgJHist[t_prime-t]))
                        Eq12[i].append(val)
                        gam = math.sqrt(2.0 * math.log(myParam.popSize))
                        Eq12_Old[i].append(math.exp((t_prime - t) * s * gam))
                        if Eq12_Old[i][len(Eq12_Old[i]) -1 ] > myParam.popSize:
                            Eq12_Old[i][len(Eq12_Old[i])-1] = 0.0
                        if Eq12[i][len(Eq12[i]) -1 ] > myParam.popSize:
                            Eq12[i][len(Eq12[i])-1] = 0.0
                avgAppear[i] += t
                found = True
                
theory_NEW_total = 0.0
for i in range(0, cellTypes):
    theoreticalAppear[i] = i * GetTau(myParam)
    theoreticalAppear_NEW[i] = theory_NEW_total + NumericalTau_TEST.SolveTauIntegral(i, myParam)
    theory_NEW_total = theoreticalAppear_NEW[i]
    avgAppear[i] /= SPD
    print("Avg Appear {0} = {1}".format(i,avgAppear[i]))
    if (i > 0):
        print("TAU {0} = {1}".format(i, avgAppear[i]-avgAppear[i-1]))
    
for i in range(2, cellTypes):
    theoreticalAppear2[i] = (i-2) * GetTau2(myParam)

found = False
target = 1.0 / (1e-7 * 100)
print("Target", target)
for i in range(0, cellTypes):
    total = 0.0
    found = False
    #if found:
     #   found = False
     #   continue
    for t in range(0, 300):
        total += Eq12[i][t]
        if total >= target and found == False:
            print("FOUND TARGET AT: " ,i, t)
            found = True

plt.figure()
plt.subplot(211)
plt.plot(avgAppear, 'o')
plt.plot(theoreticalAppear)
plt.plot(theoreticalAppear2, '--')
plt.plot(theoreticalAppear_NEW, '^')
plt.xlabel("i")
plt.ylabel("Appearance Time")

print(len(range(int(avgAppear[0]),int(avgAppear[0]) + 300)))
print(len(Eq12[0]))

plt.subplot(212)
for i in range(0, cellTypes):
   plt.plot(myHist.stepHist[0::1], myHist.histArray[i][0::1])
   plt.plot(range(int(avgAppear[i]),int(avgAppear[i]) + 300), Eq12[i], ':', linewidth = 2.0)
   plt.plot(range(int(avgAppear[i]),int(avgAppear[i]) + 300), Eq12_Old[i], '--')
plt.yscale("log")

plt.yscale("log")
plt.show()


