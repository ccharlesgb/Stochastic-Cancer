# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 12:14:12 2015

@author: Jonny
"""

import TauLeapParam
import TwoSpeciesAnalytical
import TauSolver
import anfixtime
import TwoSpecies
import wright_fisher
import TauLeap
import matplotlib.pyplot as plt
import math
import MatTools
import SimUtil


populationSize = 1e6
cellTypes = 2

#set up parameters
myParam = TauLeapParam.Params(cellTypes)
myParam.n0[0] = populationSize
myParam.SetUAll(0.05)
myParam.USE_D = False
myParam.r[0] = 1.0
myParam.r[1] = 1.0

#set up tau leaping
#Create the simulator
myTauLeaper = TauLeap.Sim(cellTypes)
myTauLeaper.params = myParam
myParam.Hook(myTauLeaper)
#Set History
tauHist = TauLeapParam.Hist(cellTypes)
#myTauLeaper.history = tauHist

myTauLeaper.n_c = 10
myTauLeaper.epsilon = 0.1
myTauLeaper.timeLimit = 10000
myTauLeaper.stopAtAppear = 0



#initialise wright fisher solver
myWF = wright_fisher.wright_fisher(cellTypes)
wfHist = TauLeapParam.Hist(cellTypes)

myWF.stopAtAppear = 0

#myWF.history = wfHist
myWF.timeLimit = 1000000
myWF.useApproxTheta = 0
myWF.params = myParam


myAnalFixtime = TwoSpeciesAnalytical.TwoSpeciesAnalytical(myParam)
#test = myAnalFixtime.getFixTime(0)
#print("Fixtime is: {0}".format(test))

DPC = 10
analytical_populations = []
min_pop = 1e1
max_pop = 1e5
for i in range(0,DPC):
    analytical_populations.append( math.ceil( SimUtil.SweepParameterLog(i,DPC,min_pop,max_pop) ))
print(analytical_populations)


wfFixtime = []
tlFixtime = []
analFixtime = []

populations = [1e1, 1e2, 1e3, 1e4, 1e5]


SPD = 5000
for i in range(0,len(populations)):
    print("Simulating - point {0} of {1}".format(i+1,len(populations)))    
    myParam.n0[0] = populations[i]
    wf_total = 0.0
    tl_total = 0.0
    for curSim in range(0,SPD):
        myTauLeaper.Simulate()
        tl_total += myTauLeaper.curTime
    
    myParam.n0[0] = populations[i]/math.sqrt(2)
    for curSime in range(0,SPD):
        myWF.Simulate()
        wf_total += myWF.curTime
    
    
    wfFixtime.append(wf_total/SPD)
    tlFixtime.append(tl_total/SPD)
    
data = dict()
data["fixation_tl"] = tlFixtime
data["analytical_populations"] = analytical_populations
data["fixation_wf"] = wfFixtime
data["populations"] = populations  



for i in range(0,len(analytical_populations)):
   print("Analytical - point {0} of {1}".format(i, len(analytical_populations)))   
   myParam.n0[0] = analytical_populations[i]
   analFixtime.append(myAnalFixtime.getFixTime(0) )
   data["fixation_analytical"] = analFixtime
   MatTools.SaveDict2(data, spd = SPD, params = myParam.GetFileString())


plt.figure()
plt.plot(populations, wfFixtime, 'o' , label = "Wright Fisher")
plt.plot(populations, tlFixtime, 'ro' , label = "Wright Fisher - No Scaling")
plt.plot(analytical_populations, analFixtime, label = "Analytical")
plt.xlabel("Population size")
plt.ylabel("Fixation Time")
plt.xscale("log")
plt.legend(loc = 2)
plt.show()
