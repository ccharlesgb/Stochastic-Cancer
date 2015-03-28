# -*- coding: utf-8 -*-
"""
Created on Mon Mar 02 15:35:16 2015

@author: Connor
"""

import matplotlib.pyplot as plt
import GillespieTauLeap_OPTON
import CarcinogenNParam_OPTON
import math
import MatTools
import wright_fisher
import SimUtil
import veusz.embed as veusz

TYPE_COUNT = 3
population = 1e6

#Create Gillespie
myParam = CarcinogenNParam_OPTON.CarcinogenNParam(TYPE_COUNT)
myGillespie = GillespieTauLeap_OPTON.Gillespie()
myGillespie.Hook(myParam) #VERY IMPORTANT
myGillespie.epsilon = 0.0001
#Create Wright Fisher
myWF = wright_fisher.wright_fisher()
myHistWF = wright_fisher.wf_hist(TYPE_COUNT)
myParamWF = wright_fisher.wright_fisher_params(TYPE_COUNT)
myWF.params = myParamWF
myWF.useApproxTheta = 0
myWF.printProgress =  1

#Set sim params
myGillespie.timeLimit = 40000
myGillespie.RECORD_TAU_INFO = 0
myGillespie.printProgress = 1
myGillespie.stopAtAppear = 1

#Set Parameters
myParam.n0[0] = population
myParam.d = 100
#Set mutation and fitness
s = 0.0
for i in range(0,TYPE_COUNT):
    myParam.r[i] = math.pow(1.0 + s, i)
    myParam.u[i] = 1e-5

def SetWFParams(gilParam):
    #Set Sim params
    myWF.stepLimit = myGillespie.timeLimit
    
    myParamWF.d = gilParam.d;
    myParamWF.iN[0] = gilParam.n0[0]
    for i in range(0, TYPE_COUNT):
        myParamWF.r[i] = gilParam.r[i]
        myParamWF.u[i] = gilParam.u[i]

SetWFParams(myParam)

SDP = 25
pointCount = 8

dataX = []
dataY_Gill = []
dataY_Tau = []
dataY_WF = []
dataY_err = []

minN = 1e4
maxN = 1e1

minU = 1e-4
maxU = 1e-2

minR = 0.3
maxR = 3.0

for i in range(0, pointCount):
    myParam.n0[0] = int(SimUtil.SweepParameterLog(i,pointCount, minN, maxN))
    print("N = {0}".format(myParam.n0[0]))
    dataX.append(myParam.n0[0])
    #myParam.r[1] = SimUtil.SweepParameter(i, pointCount, minR, maxR)
    #print("R = {0}".format(myParam.r[1]))
    #dataX.append(myParam.r[1])
    
    SetWFParams(myParam)

    print("WF "),
    wfRes = myWF.SimulateBatch(SDP)
    print("GI "),
    gillRes = myGillespie.SimulateBatch(SDP)
    print("TAU"),
    myGillespie.epsilon = 0.1
    tauRes = myGillespie.SimulateBatch(SDP)
    
    dataY_Gill.append(gillRes.avgFixTime)
    dataY_Tau.append(tauRes.avgFixTime)
    dataY_WF.append(wfRes.avgFixTime)
    dataY_err.append(float(wfRes.avgFixTime - gillRes.avgFixTime))
    
plt.subplot(211)
plt.plot(dataX, dataY_Gill, '^-')
plt.plot(dataX, dataY_Tau, '+-')
plt.plot(dataX, dataY_WF ,'o-')
plt.xscale("log")

plt.subplot(212)
plt.plot(dataX, dataY_err)
plt.plot()
plt.xscale("log")

file_name = "MethodComparison_SweepN_DPC={0}_SDP_{1}_CT_{2}_U_{3}_N_{4}_{5}".format(pointCount, SDP, TYPE_COUNT, myParam.u[0], minN, maxN)

data = dict()
data["N"] = dataX
data["Nt_{0}_GI".format(TYPE_COUNT)]  = dataY_Gill
data["Nt_{0}_TA".format(TYPE_COUNT)] = dataY_Tau
data["Nt_{0}_WF".format(TYPE_COUNT)] = dataY_WF

MatTools.SaveDict(file_name, data)

g = veusz.Embedded('window title')
g.EnableToolbar()

g.To( g.Add('page') )
g.To( g.Add('graph') )
g.Add('xy', marker='tiehorz', MarkerFill__color='green')
g.SetData('x', dataX)
g.SetData('y', dataY_Gill)

