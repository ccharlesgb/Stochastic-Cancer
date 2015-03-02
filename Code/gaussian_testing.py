# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 15:00:22 2015

@author: Jonny
"""

import wright_fisher
import matplotlib.pyplot as plt
import math
import gaussian_fitter

#set up simulation
cellTypes = 21
population = 1e9

myWF = wright_fisher.wright_fisher(cellTypes)
myHist = wright_fisher.wf_hist(cellTypes)

myWF.u = 1e-9
#myWF.u = 1e-2 / population

myWF.iN[0] = population
myWF.history = myHist
myWF.popSize = population

myWF.stepLimit = 10000
myWF.useApproxTheta = 0
myWF.d = 200

s = 0.01
adv = 0
for i in range(0,cellTypes):
    if i % 2 == 0:
        myWF.r[i] = math.pow(1.0 + s, adv)
        adv += 1
    else:
        myWF.r[i] = math.pow(1.0 + s, adv)
        adv += 1


myWF.Simulate()

#set up fitter
myFitter = gaussian_fitter.gaussian_fitter()
myFitter.return_data = 1
plt.figure()
step = myWF.curStep / 6
first = int(step)
end = int(myWF.curStep - step) 

for t in range(first, end, step):
    dataX = []
    dataY = []
    for i in range(0, cellTypes):
        dataX.append(i)
        dataY.append(myHist.histArray[i][t])     
    
    myFitter.SetData(dataX,dataY)
    myFitter.FitData()
    #print("The standard dev at time {0} is: {1}".format(t, myFitter.sdev))
    
    plt.plot(dataX,dataY, 'o-',label = "Simulation")
    #print("dataX is {0}".format(dataX))
    #print("dataY is {0}".format(myFitter.fitted_y_data))
    plt.plot(dataX, myFitter.fitted_y_data, 'x-',label = "Fitted")
    
plt.yscale("log")
plt.xlabel("Number of Mutations")
plt.ylabel("Cell count")
plt.show()
