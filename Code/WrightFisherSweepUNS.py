# -*- coding: utf-8 -*-
"""
Created on Thu Feb 26 13:32:23 2015
@author: Connor
"""

import wright_fisher
import matplotlib.pyplot as plt
import math
import SimUtil
import time
import MatTools
import TauSolver
import TauLeapParam

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

myParam.SetUAll(1e-7)

s = 0.01
for i in range(0,cellTypes):
    myParam.r[i] = math.pow(1.0 + s, i)

mySolver = TauSolver.Solver(myParam)    

minN = 1e6
maxN = 1e9

minU = 1e-8
maxU = 1e-5

minS = 1e-4
maxS = 1e-1

SDP = 1
PointCount = 4
PointCountTheory = 10

NdataX = []
NdataY = []
NdataX_orig = []
NdataY_orig = []

UdataX = []
UdataY = []
UdataX_orig = []
UdataY_orig = []

SdataX = []
SdataY = []
SdataX_orig = []
SdataY_orig = []


for p in range(0,PointCount):
    startTime = time.clock()
    print("Current Data Point = {0}/{1} ({2}%)".format(p + 1, PointCount, 100.0 * float(p)/(PointCount-1)))
    #Sweep N
    myParam.n0[0] = int(SimUtil.SweepParameterLog(p,PointCount, minN, maxN))
    myParam.SetUAll(1e-7)
    myParam.SetCompoundFitness(1e-2)
    print("N", myParam.n0[0])
    res = myWF.SimulateBatch(SDP)
    NdataX.append(myParam.n0[0])
    NdataY.append(res.avgFixTime)
    
    #Sweep U
    myParam.n0[0] = 1e9
    myParam.SetCompoundFitness(1e-2)
    myParam.SetUAll(SimUtil.SweepParameterLog(p,PointCount, minU, maxU))
    myParam.Reset()
    res = myWF.SimulateBatch(SDP)
    UdataX.append(myParam.u[0])
    UdataY.append(res.avgFixTime)
    
    #Sweep S
    myParam.n0[0] = 1e9
    myParam.SetUAll(1e-7)
    s = SimUtil.SweepParameterLog(p,PointCount, minS, maxS)
    myParam.SetCompoundFitness(s)
    myParam.Reset()
    res = myWF.SimulateBatch(SDP)
    SdataX.append(s)
    SdataY.append(res.avgFixTime)
    print("Complete (Took {:.{s}f} seconds)".format(time.clock() - startTime, s=1))    
    

for p in range(0, PointCountTheory):
    #Sweep N
    myParam.n0[0] = int(SimUtil.SweepParameterLog(p,PointCountTheory, minN, maxN))
    myParam.SetUAll(1e-7)
    myParam.SetCompoundFitness(1e-2)
    myParam.Reset()
    theory1 = mySolver.GetWaitingTimeOriginal(cellTypes - 1)
    NdataX_orig.append(myParam.n0[0])
    NdataY_orig.append(theory1)
    #Sweep U
    myParam.n0[0] = 1e9
    myParam.SetCompoundFitness(1e-2)
    myParam.SetUAll(SimUtil.SweepParameterLog(p,PointCountTheory, minU, maxU))
    myParam.Reset()
    theory1 = mySolver.GetWaitingTimeOriginal(cellTypes - 1)
    UdataX_orig.append(myParam.u[0])
    UdataY_orig.append(theory1)
    #Sweep S
    myParam.n0[0] = 1e9
    myParam.SetUAll(1e-7)
    s = SimUtil.SweepParameterLog(p,PointCountTheory, minS, maxS)
    myParam.SetCompoundFitness(s)
    myParam.Reset()
    theory1 = mySolver.GetWaitingTimeOriginal(cellTypes - 1)
    SdataX_orig.append(s)
    SdataY_orig.append(theory1)


plt.figure()
plt.subplot(131)
plt.plot(NdataX,NdataY, 'o')
plt.plot(NdataX_orig,NdataY_orig)
plt.xscale("log")
plt.xlabel("N")
plt.ylabel("t_{0}".format(cellTypes - 1))
plt.xlim(minN, maxN)
plt.subplot(132)
plt.plot(UdataX,UdataY, 'o')
plt.plot(UdataX_orig,UdataY_orig)
plt.xscale("log")
plt.xlabel("U")
plt.ylabel("t_{0}".format(cellTypes - 1))
plt.xlim(minU, maxU)
plt.subplot(133)
plt.plot(SdataX,SdataY, 'o')
plt.plot(SdataX_orig,SdataY_orig)
plt.xscale("log")
plt.xlabel("S")
plt.ylabel("t_{0}".format(cellTypes - 1))
plt.xlim(minS, maxS)



plt.show()

data = dict()
data["NX"] = NdataX
data["Nt_20"] = NdataY
data["NX_orig"] = NdataX_orig
data["Nt_20_orig"] = NdataY_orig
data["UX"] = UdataX
data["Ut_20"] = UdataY
data["UX_orig"] = UdataX_orig
data["Ut_20_orig"] = UdataY_orig
data["SX"] = SdataX
data["St_20"] = SdataY
data["SX_orig"] = SdataX_orig
data["St_20_orig"] = SdataY_orig


MatTools.SaveDict2(data,N='N',SDP = SDP, DPC = PointCount, PARAM = myParam.GetFileString())
