# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 18:53:15 2015

@author: Jonny
"""

import wright_fisher
import matplotlib.pyplot as plt
import SimUtil
import math

def GetT_K(param):
    summation = 0.0
    s = param.r[1] - param.r[0]
    for i in range(0, param.cellTypes):
        top = math.pow(math.log(s / (param.u[i] * param.d)), 2.0)
        bottom = 2.0 * s * math.log(param.popSize)
        summation += (top / bottom)
    return summation;

import wright_fisher
import matplotlib.pyplot as plt
import math
import SimUtil
import NumericalTau_TEST

def GetTau2(param):
    s = param.r[1] - param.r[0]
    logs = math.log(s / (param.u[0] * param.d)*math.sqrt(2.0*math.log(param.popSize)))
    top = math.pow(logs,2.0)
    bottom = 2.0 * s * math.log(param.popSize)
    return top/bottom

def GetTau(param):
    s = param.r[1] - param.r[0]
    top = math.pow(math.log(s / (param.u[0] * param.d)),2.0)
    bottom = 2.0 * s * math.log(param.popSize)
    return top/bottom

#Get the number of cells expected from one step.
def SolveXJ0(j, param):
    #print("GETXJ0", j)
    if j == 0:
        return 1.0
    if j == 1:
        return 1e-5
    fac = param.u[0] * param.d
    t = 0.0
    maxT = 1.0
    step = 0.1
    
    result = 0.0
    while (t < maxT):
        result = result + step * GetXJ(t,j-1,param)
        t = t + step
    #print("X_{0}(0) = {1}".format(j, result * fac))
    return max(result * fac, 1/param.popSize)

def GetXJ(t,j,param):
    #print("GETXJ", j)
    xj = 0.0    
    s = param.r[1] - param.r[0]
    if j == 0:
        return 1.0
    x_j0 = SolveXJ0(j, param)
    x_j_1 = SolveXJ0(j-1,param)
    #((u * d * x_j_1 + (s * (i-jbar) + u*d) * start/1e9) * math.exp((s-x_j*u*d) * (i - jbar) * T) - u * d * x_j_1) / (s * (i - jbar) + u*d)
    #xj = ((param.u[0] * param.d * x_j_1 + (s * j + u*d) * x_j0) * math.exp((s-x_j*u*d) * (i - jbar) * T) - u * d * x_j_1) / (s * j + u*d)
    xj = 1.0/(s * j) * ((param.u[0] * param.d * x_j_1 + s * j * x_j0) * math.exp(s * j * t) - param.u[0] * param.d * x_j_1)
    #print("X {0}(t={1}) IS {1}".format(j,t, xj))
    return xj
    
def GetTau3(j,param):
    fac = param.popSize * param.u[0] * param.d
    t = 0.0
    maxT = 1.0
    step = 1e-3
    
    result = 0.0
    while (result < 1.0):
        result += step * fac * GetXJ(t,j,param)
        if result >= 1.0:
            break;
        t = t + step

    return t

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

s = 0.01
for i in range(0,cellTypes):
    myParam.r[i] = math.pow(1.0 + s, i)

uPowMin = 1.8
uPowMax = 2.2

SPD = 1
DPC = 5

s = 0.01
for i in range(0,cellTypes):
    myParam.r[i] = math.pow(1.0 + s, i)

myWF.reset()

SPD = 10
DPC = 5

dataX = []
dataY = []
dataY_Theory = []

s = 0.01
for i in range(0,cellTypes):
    myParam.r[i] = math.pow(1.0 + s, i)

for curPoint in range(0,DPC):
    uPow = SimUtil.SweepParameter(curPoint, DPC, uPowMin, uPowMax)
    for i in range(0,cellTypes):
        myParam.u[i] = 1e-10 * math.pow(uPow,i)
    
    print("U: ", myParam.u[0], myParam.u[cellTypes - 1])    

    dataX.append(uPow)

    res = myWF.SimulateBatch(SPD)
    dataY.append(res.avgFixTime)
    dataY_Theory.append(GetT_K(myParam))
    
plt.plot(dataX, dataY, 'o')
plt.plot(dataX, dataY_Theory)
plt.xlabel("uPow")
plt.ylabel("Steps to fixation")

