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
import TauLeapParam
import matplotlib.colors as colors
import matplotlib.cm as cm
        
errormap = colors.ListedColormap([[0.1, 0.2, 1.0],
                       [0.4, 0.8, 1.0],
                       [0.5, 1.0, 0.3],
                       [1.0, 0.5, 0.1],
                       [1.0, 0.18, 0.15]])

cm.register_cmap('Errors', cmap = errormap)


cellTypes = 21
population = 1e9

myWF = wright_fisher.wright_fisher(cellTypes)
myParam = TauLeapParam.Params(cellTypes)

myWF.stopAtAppear = 1
myWF.timeLimit = 100000000
myWF.useApproxTheta = 0
myWF.params = myParam

myParam.d = 100
myParam.n0[0] = population
myParam.SetUAll(1e-7)

s = 0.01
for i in range(0,cellTypes):
    myParam.r[i] = math.pow(1.0 + s, i)

mySolver = TauSolver.Solver(myParam)

mapSize = 8

minS = 1e-4
maxS = 1e-1

minU = 1e-7
maxU = 1e-4

SDP = 30

sHist = []
uHist = []

SU_CUTOFF = 0 #Set to 0 to disable cutoff
SU_ROUGH = 2 #Cut the SDP here as it takes ages
SDP_ROUGH = 5
DONT_SIM = False

#Set up our arrays
avgFixTime = np.zeros((mapSize, mapSize))
fixTime1 = np.zeros((mapSize, mapSize))
fixTime2 = np.zeros((mapSize, mapSize))
fixTime3 = np.zeros((mapSize, mapSize))
fixTime4 = np.zeros((mapSize, mapSize))
errFixTime1 = np.zeros((mapSize, mapSize))
errFixTime2 = np.zeros((mapSize, mapSize))
errFixTime3 = np.zeros((mapSize, mapSize))
errFixTime4 = np.zeros((mapSize, mapSize))

#Create the map
for iS in range(0, mapSize):
    for iU in range(0, mapSize):
        s = SimUtil.SweepParameterLog(iS, mapSize, minS,maxS)
        u = SimUtil.SweepParameterLog(iU, mapSize, minU,maxU)
        sHist.append(math.log(s, 10.0))
        uHist.append(math.log(u, 10.0))
        print("s = {0} u = {1}".format(s,u))        
        
        myParam.SetUAll(u)
        for i in range(0,cellTypes):
            myParam.r[i] = math.pow(1.0 + s, i)
        mySolver.CacheX0()
        
        if DONT_SIM == True  or (iS < SU_CUTOFF and iU < SU_CUTOFF):
            res = myWF.SimulateBatch(0)
            res.avgFixTime = 1e9
        elif (iS < SU_ROUGH and iU < SU_ROUGH): #Faster simulation here
            res = myWF.SimulateBatch(SDP_ROUGH)
        else:
            res = myWF.SimulateBatch(SDP)
        #res.avgFixTime = random.randrange(800, 5000)
        avgFixTime[iU, iS] = res.avgFixTime
        
        theory1 = mySolver.GetWaitingTimeOriginal(cellTypes - 1)
        theory2 = mySolver.GetWaitingTimeNeglect(cellTypes - 1)
        theory3 = mySolver.GetWaitingTimeModelNew(cellTypes - 1)
        theory4 = mySolver.GetWaitingTimeModel(cellTypes - 1)
        
        fixTime1[iS,iU] = theory1
        fixTime2[iS,iU] = theory2
        fixTime3[iS,iU] = theory3
        fixTime4[iS,iU] = theory4
        
        errFixTime1[iS, iU] = (theory1 - res.avgFixTime) / res.avgFixTime
        errFixTime2[iS, iU] = (theory2 - res.avgFixTime) / res.avgFixTime 
        errFixTime3[iS, iU] = (theory3 - res.avgFixTime) / res.avgFixTime   
        errFixTime4[iS, iU] = (theory4 - res.avgFixTime) / res.avgFixTime   
        

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
plt.pcolormesh(x,y,errFixTime1, cmap = errormap)
plt.colorbar()
plt.xlabel("S")
plt.ylabel("U")
plt.clim(-1.0,1.0)
plt.title("Original")

plt.subplot(222)
plt.pcolormesh(x,y,errFixTime2, cmap = errormap)
plt.colorbar()
plt.xlabel("S")
plt.ylabel("U")
plt.clim(-1.0,1.0)
plt.title("Neglect")

plt.subplot(223)
plt.pcolormesh(x,y,errFixTime3, cmap = errormap)
plt.colorbar()
plt.xlabel("S")
plt.ylabel("U")
plt.clim(-1.0,1.0)
plt.title("Model Iterative")

plt.subplot(224)
plt.pcolormesh(x,y,errFixTime4, cmap = errormap)
plt.colorbar()
plt.xlabel("S")
plt.ylabel("U")
plt.clim(-1.0,1.0)
plt.title("Mutational Correction")

data1 = dict()
data1["xCoords"] = x
data1["yCoords"] = y
data1["time"] = avgFixTime

data2 = dict()
data2["xCoords"] = x
data2["yCoords"] = y
data2["fixTimeOrig"] = fixTime1
data2["fixTimeNegl"] = fixTime2
data2["fixTimeModN"] = fixTime3
data2["fixTimeModA"] = fixTime4
data2["errFixTimeOrig"] = errFixTime1
data2["errFixTimeNegl"] = errFixTime2
data2["errFixTimeModN"] = errFixTime3
data2["errFixTimeModA"] = errFixTime4
    
if DONT_SIM == False:
    MatTools.SaveDict2(data1, TIMES = "", SDP = SDP, SIZE = mapSize, N = population)
MatTools.SaveDict2(data2, PREDICT = "", SDP = SDP, SIZE = mapSize, N = population, DONT_SIM = DONT_SIM)
