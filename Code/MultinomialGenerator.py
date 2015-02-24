# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 14:52:31 2015

@author: Connor
"""

import math
import random
import matplotlib.pyplot as plt
import numpy as np
import MatTools

pVec = []
cumPVec = []

CELL_COUNT = int(1e9)
X = []

TYPE_COUNT = 21

total = 0.0
for i in range(0, TYPE_COUNT):
    pVec.append(0.0)
    total += pVec[i]
    X.append(0)

STEPS = 30000

s = 1e-2
u = 1e-7

X[0] = CELL_COUNT

r = []

sigma = 0.01

for i in range(0,TYPE_COUNT):
    r.append(random.gauss(math.pow(1.0 + 0.01, i), sigma))
    #r.append(math.pow(1.0 + 0.01, i))

print(r)

tHist = []
histArray = dict()
for i in range(0, TYPE_COUNT):
    histArray[i] = []

step = 0
fixed = 0
while (step < STEPS and fixed != 1):
    #print(step)
    total = 0.0
    # getting theta_i    
    for i in range(0, TYPE_COUNT):
        avgFit = 0.0
        for l in range(0, TYPE_COUNT):
            avgFit += r[l]*X[l]
        #print(avgFit)
        theta_i = (r[i] * X[i])/avgFit
        if i > 0:
            theta_i += float(u * (TYPE_COUNT - i + 1) * r[i-1]*X[i-1] )/ avgFit
            
        pVec[i] = theta_i
        total += pVec[i]


    # get proba vector
    for i in range(0, TYPE_COUNT):
        pVec[i] /= total
        if pVec[i] == 1.0:
            print("FIXATION")
            fixed = 1

    X = np.random.multinomial(CELL_COUNT, pVec)
    
    if X[TYPE_COUNT - 1] >= 1:
        print("GOT MUTANT")
        fixed = 1
    
    tHist.append(step)
    for i in range(0, TYPE_COUNT):
        histArray[i].append(X[i])
        
    step += 1
    '''
data = dict()    
    
for t in range(0, step, step / 5):
    dataX = []
    dataY = []
    for i in range(0, TYPE_COUNT):
        dataX.append(i)
        dataY.append(histArray[i][t])
     
    data["X_{0}".format(t)] = dataX
    data["Y_{0}".format(t)] = dataY
    plt.plot(dataX,dataY, '-o')
    plt.yscale("log")
        '''
        
for i in range(0, TYPE_COUNT):
    plt.plot(tHist, histArray[i])
    plt.yscale("log")

file_name = "GaussianDistributions"
'''
MatTools.SaveDict(file_name, data)'''