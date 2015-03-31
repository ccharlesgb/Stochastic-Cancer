# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 15:00:22 2015

@author: Jonny
"""

import wright_fisher
import matplotlib.pyplot as plt
import math
import gaussian_fitter
import TauLeapParam

#set up simulation
cellTypes = 21
population = 1e9

myWF = wright_fisher.wright_fisher(cellTypes)
myWF_params = TauLeapParam.Params(cellTypes)
myWF.params = myWF_params

myHist = TauLeapParam.Hist(cellTypes)

myWF_params.u = [1e-9]*cellTypes
#myWF.u = 1e-2 / population

myWF_params.n0[0] = population
myWF.history = myHist
myWF.popSize = population

myWF.timeLimit = 1e4
myWF.useApproxTheta = 0
myWF.d = 200

s = 0.01
for i in range(0,cellTypes):
    if i % 2 == 0:
        myWF_params.r[i] = math.pow(1.0 + s, i)
    else:
        myWF_params.r[i] = math.pow(1.0 + s, i)

#set up fitter
myFitter = gaussian_fitter.gaussian_fitter()
#myFitter.return_data = 1



SPD = 100
DPC = 6
fits_per_sim = 3

sigma = []
mut_rate = [1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4]
min_u = 1e-10
max_u = 1e-5


step = myWF.curTime / fits_per_sim
first = int(step)
end = int(myWF.curTime) 

for curPoint in range(0,len(mut_rate)):
    print("Datapoint {0} of {1}".format(curPoint + 1, DPC))        
    for mut in range(0,cellTypes):    
        myWF_params.u[mut] = mut_rate[curPoint]
    sigma_term = 0.0    
    errorcount = 0
    for sim in range(0,SPD):
        myWF.Simulate()
        t = int(myWF.curTime / 2)
        '''        
        for t in range(int(myWF.curStep / fits_per_sim), int(myWF.curStep), int(myWF.curStep / fits_per_sim)):
            dataX = []
            dataY = []
            for i in range(0, cellTypes):
                dataX.append(i)
                dataY.append(myHist.histArray[i][t])     
        ''' 
        dataX = []
        dataY = []        
        for i in range(0, cellTypes):
                dataX.append(i)
                dataY.append(myHist.histArray[i][t])     
        
        myFitter.SetData(dataX, dataY)
        myFitter.FitData()
        sigma_term += abs(myFitter.sdev)
        errorcount += myFitter.errorcount
        '''
        plt.plot(dataX,dataY, 'o-',label = "Simulation")
        print("dataX is {0}".format(dataX))
        print("dataY is {0}".format(myFitter.fitted_y_data))
        plt.plot(dataX, myFitter.fitted_y_data, 'x-',label = "Fitted")
        '''
    sigma.append(sigma_term/(SPD-errorcount))

plt.plot(mut_rate, sigma)
plt.xlabel("Mutational Rate")
plt.ylabel("Fitted Gaussian Sigma")
plt.xscale("log")
plt.show()

'''        
plt.yscale("log")
plt.xlabel("Number of Mutations")
plt.ylabel("Cell count")
plt.show()
'''