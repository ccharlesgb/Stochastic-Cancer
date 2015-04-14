# -*- coding: utf-8 -*-
"""
Created on Wed Apr 08 15:29:41 2015

@author: Connor
"""

'''
IDEA:
Set mutation rate to u_i = initial mutation rate
Model effect of a gene that causes chromosonal instability mutating as a step
function. Ie once this gene has mutated it causes u to go from u_i to u_f which
would be much higher.
Model at which cell type causes the worst increase in fixation time and by how
much.
'''
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

SDP = 20
PointCount = cellTypes - 1

dataX = []
dataY = []
dataY_anal1 = []

u_before = 1e-7
u_after = 1e-5

for p in range(0,PointCount):
    startTime = time.clock()
    print("Current Data Point = {0}/{1} ({2}%)".format(p + 1, PointCount, 100.0 * float(p)/(PointCount-1)))
    
    for i in range(1, cellTypes):
        if (p + 1) > i:
            myParam.SetU(i, u_before)
        else:
            myParam.SetU(i, u_after)
    print(myParam.u)

    res = myWF.SimulateBatch(SDP)
    dataX.append(p)
    dataY.append(res.avgFixTime)
    
    mySolver.CacheX0()    
    
    theory4 = mySolver.GetWaitingTimeModel(cellTypes - 1)
    
    dataY_anal1.append(theory4)
    
    print("Complete (Took {:.{s}f} seconds)".format(time.clock() - startTime, s=1))    

plt.figure()
plt.plot(dataX,dataY, 'o')
plt.plot(dataX,dataY_anal1, ':')
plt.xlabel("gene_on")
plt.ylabel("t_{0}".format(cellTypes - 1))

plt.show()

data = dict()
data["N"] = dataX
data["Nt_20"] = dataY
data["Nt_20_anal1"] = dataY_anal1

MatTools.SaveDict2(data,SDP = SDP, DPC = PointCount, PARAM = myParam.GetFileString())
