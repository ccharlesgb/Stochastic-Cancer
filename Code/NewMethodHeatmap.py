# -*- coding: utf-8 -*-

import numpy as np

import wright_fisher
import matplotlib.pyplot as plt
import math
import SimUtil
import time
import MatTools
import TauSolver
import random


cellTypes = 21
population = 1e9

myWF = wright_fisher.wright_fisher()
myParam = wright_fisher.wright_fisher_params(cellTypes)

myWF.stopAtAppear = 1
myWF.stepLimit = 10000000
myWF.useApproxTheta = 0
myWF.params = myParam

myParam.d = 100
myParam.iN[0] = population
myParam.u[0] = 1e-7

s = 0.01
for i in range(0,cellTypes):
    myParam.r[i] = math.pow(1.0 + s, i)

mySolver = TauSolver.Solver(myParam)

mapSize = 8

minS = 1e-4
maxS = 1e-1

minU = 1e-8
maxU = 1e-5

SDP = 30

avgFixTime = np.zeros((mapSize, mapSize))

errFixTime1 = np.zeros((mapSize, mapSize))
errFixTime2 = np.zeros((mapSize, mapSize))
errFixTime3 = np.zeros((mapSize, mapSize))

#xticks = np.arange(0, mapSize, 5.0/(mapSize-1.0))
#xlabels = np.arange(minS, maxS, 1.0/(mapSize-1.0))

#yticks = np.arange(0, mapSize, 5.0/(mapSize-1.0))
#ylabels = np.arange(minU, maxU, 1.0/(mapSize-1.0))

sHist = []
uHist = []

for iS in range(0, mapSize):
    for iU in range(0, mapSize):
        s = SimUtil.SweepParameterLog(iS, mapSize, minS,maxS)
        u = SimUtil.SweepParameterLog(iU, mapSize, minU,maxU)
        sHist.append(math.log(s, 10.0))
        uHist.append(math.log(u, 10.0))
        print("s = {0} u = {1}".format(s,u))        
        
        myParam.u[0] = u
        for i in range(0,cellTypes):
            myParam.r[i] = math.pow(1.0 + s, i)
        mySolver.CacheX0()
        
        res = myWF.SimulateBatch(SDP)
        #res.avgFixTime = random.randrange(800, 5000)
        avgFixTime[iU, iS] = res.avgFixTime
        
        theory1 = mySolver.GetWaitingTimeOriginal(cellTypes)
        theory2 = mySolver.GetWaitingTimeNeglect(cellTypes)
        theory3 = mySolver.GetWaitingTime(cellTypes)
        
        errFixTime1[iS, iU] = (theory1 - res.avgFixTime) / res.avgFixTime
        errFixTime2[iS, iU] = (theory2 - res.avgFixTime) / res.avgFixTime 
        errFixTime3[iS, iU] = (theory3 - res.avgFixTime) / res.avgFixTime   
        

minSHist = round(min(sHist))
maxSHist = round(max(sHist))

minUHist = round(min(uHist))
maxUHist = round(max(uHist))
print(minUHist, maxUHist)
print(minSHist, maxSHist)

dy = abs(maxUHist - minUHist) / mapSize
dx = abs(maxSHist - minSHist) / mapSize

print("MAXU", maxUHist)
print("MAXS", maxSHist)

print("DX IS", dx)
y, x = np.mgrid[slice(minUHist, maxUHist + dy, dy), slice(minSHist, maxSHist + dx, dx)]

#print(y)
#print(x)

plt.figure()
plt.subplot(221)
plt.pcolormesh(x,y,avgFixTime)
plt.colorbar()
plt.xlabel("S")
plt.ylabel("U")

plt.subplot(222)
plt.pcolormesh(x,y,errFixTime1)
plt.colorbar()
plt.xlabel("S")
plt.ylabel("U")
plt.clim(-1.0,1.0)

plt.subplot(223)
plt.pcolormesh(x,y,errFixTime2)
plt.colorbar()
plt.xlabel("S")
plt.ylabel("U")
plt.clim(-1.0,1.0)

plt.subplot(224)
plt.pcolormesh(x,y,errFixTime3)
plt.colorbar()
plt.xlabel("S")
plt.ylabel("U")
plt.clim(-1.0,1.0)

file_name = "NewMethodErrorMap_SDP={0}_Size={1}_N={2}".format(SDP, mapSize, population)

data = dict()
data["xCoords"] = x
data["yCoords"] = y
data["time"] = avgFixTime
data["errFixTime1"] = errFixTime1
data["errFixTime2"] = errFixTime2
data["errFixTime3"] = errFixTime3
    
MatTools.SaveDict(file_name, data)