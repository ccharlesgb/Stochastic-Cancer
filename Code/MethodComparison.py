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

TYPE_COUNT = 11
population = 1e6

myParam = CarcinogenNParam_OPTON.CarcinogenNParam(TYPE_COUNT)
#Create the simulator
myGillespie = GillespieTauLeap_OPTON.Gillespie()
myGillespie.Hook(myParam) #VERY IMPORTANT
#Set History
myHist = CarcinogenNParam_OPTON.CarcinogenNHist(TYPE_COUNT)
#myGillespie.history = myHist

myGillespie.timeLimit = 40000
myGillespie.RECORD_TAU_INFO = 0
myGillespie.printProgress = 1
myGillespie.stopAtAppear = 1
myParam.n0[0] = population

myParam.d = 100

myWF = wright_fisher.wright_fisher()
myHistWF = wright_fisher.wf_hist(TYPE_COUNT)
myParamWF = wright_fisher.wright_fisher_params(TYPE_COUNT)

myParamWF.d = 100

myParamWF.iN[0] = int(population) / math.sqrt(2.0)
#Do this because WF can spawn itself but Tau cannot

#myWF.history = myHistWF
myWF.stepLimit = 40000
myWF.useApproxTheta = 0
myWF.params = myParamWF
myWF.printProgress =  1

#Set mutation and fitness
s = 0.01
for i in range(0,TYPE_COUNT):
    myParam.r[i] = math.pow(1.0 + s, i)
    myParam.u[i] = 1e-3
    myParamWF.r[i] = math.pow(1.0 + s, i)
    myParamWF.u[i] = 1e-3
    
myGillespie.epsilon = 0.05

SDP = 2500
pointCount = 8

dataX = []
dataY_Gill = []
dataY_WF = []

minN = 1e4
maxN = 1e1

for i in range(0, pointCount):
    myParam.n0[0] = int(SimUtil.SweepParameterLog(i,pointCount, minN, maxN))
    myParamWF.iN[0] = int(myParam.n0[0] / math.sqrt(2.0))
    dataX.append(myParam.n0[0])
    print("N = {0}".format(myParam.n0[0]))
    print("WF"),
    wfRes = myWF.SimulateBatch(SDP)
    print("GI"),
    gillRes = myGillespie.SimulateBatch(SDP)
    
    dataY_Gill.append(gillRes.avgFixTime)
    dataY_WF.append(wfRes.avgFixTime)
    
plt.plot(dataX, dataY_Gill)
plt.plot(dataX, dataY_WF)
plt.xscale("log")

file_name = "MethodComparison_SweepN_DPC={0}_SDP_{1}_CT_{2}_U_{3}_N_{4}_{5}".format(pointCount, SDP, TYPE_COUNT, myParam.u[0], minN, maxN)

data = dict()
data["N"] = dataX
data["Nt_{0}_GI".format(TYPE_COUNT)]  = dataY_Gill
data["Nt_{0}_WF".format(TYPE_COUNT)] = dataY_WF

MatTools.SaveDict(file_name, data)