# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 18:41:05 2015

@author: Connor
"""

import math
import matplotlib.pyplot as plt
import wright_fisher

cellTypes = 21
population = 1e9

myParam = wright_fisher.wright_fisher_params(cellTypes)

myParam.d = 100

myParam.iN[0] = population

myParam.u[0] = 1e-7

s = 0.01
for i in range(0,cellTypes):
    myParam.r[i] = math.pow(1.0 + s, i)

myParam.Reset()

def TauNought(j, param): 
    return max(param.u[0] * param.d * TauInt(j-1, 1.0, param), 1.0/param.popSize)

def TauInt(j, tau, param):
    s = param.r[1] - param.r[0]
    u = param.u[0]
    d = param.d
    N = param.popSize        
    u = 0.0
    if j < 0:
        return 0.0
    if j == 0:
        return tau
    x_j_0 = TauNought(j, param)
    x_j_1 = TauNought(j - 1, param)   
    
    #print("1x_{0}_0 = {1}".format(j,x_j_0))
    #print("2x_{0}_0 = {1}".format(j-1,x_j_1))
    '''if j == 0:
        x_j_0 = 1.0
        x_j_1 = 0.0
    if j == 1:
        x_j_0 = 1e4 / N
        x_j_1 = 1.0    
    if j == 2:
        x_j_1 = 1e4 / N '''

    gamma = math.sqrt(2.0 / (s * tau) * math.log(1.0 / (x_j_1)))
    #a = s * gamma - u * d
    b = u * d * x_j_1
    c = x_j_0
    
    a = s * gamma
    
    #top = (a*c + b)*math.exp(a * tau) - a * b * tau - a *c - b
    top = (a*c + b)*math.exp(a * tau) - a * b * tau - a * c
    top = top / (a*a)
    top = top
    
    return top

def SolveTauIntegral(j, param):
    foundSol = False
    epsilon = 1e-14 #Min value f prime
    tolerance = 1e-5 # accuracy required
    maxIter = 200
    
    delta = 0.01
    x0 = 10.0
    for i in range(0, maxIter):
        y = TauInt(j, x0, param) - 1.0 / (param.popSize * param.u[0] * param.d)
        y_prime = ((TauInt(j, x0 + delta, param) - 1.0 / (param.popSize * param.u[0] * param.d)) - y) / delta
        #print("x0 = {0} y' = {1}".format(x0, y_prime))
        if abs(y_prime) < epsilon:
            break #Denominator too small
        x1 = min(x0 - y / y_prime, 1000.0)
        
        if abs(x1 - x0)/abs(x1) < tolerance:
            foundSol = True
            break
        x0 = x1
    if foundSol == True:
        return x0
    else:
        return -1.0
            
dataX = []
dataY = []

PC = 15
for i in range(0, PC):
    dataX.append(i / 1000.0)
    dataY.append(TauInt(0, dataX[i],myParam))
    

#plt.plot(dataX,dataY)

#print("ANSWER ", SolveTauIntegral(0, myParam))
#print("ANSWER ", SolveTauIntegral(1, myParam))
#print("ANSWER ", SolveTauIntegral(2, myParam))

def GetXJ(t,j, param):
    s = param.r[1] - param.r[0]
    u = param.u[0]
    d = param.d
    N = param.popSize    
    
    x_j_0 = 1.0 / N
    x_j_1 = 1.0 / N     
    if j == 0:
        x_j_0 = 1.0 - 1.0/N
        x_j_1 = 0.0
    if j == 1:
        x_j_0 = 1e4 / N
        x_j_1 = 1.0 - 1.0/N
    if j == 2:
        x_j_1 = 1e4 / N 
    t = float(t)
    tau_epsilon = 1e-10
    gamma = math.sqrt((2.0 / (s * (t+tau_epsilon))) * math.log(1.0 / (x_j_1)))
    #gamma = math.sqrt(2 * math.log(1.0 / (x_j_0)))
    a = s * gamma
    b = u * d * x_j_1
    c = x_j_0
    a = s * gamma
    
    res = (b + c*a)*math.exp(a*t) - b
    return 1.0/a * res
    
    
    
    
    
    