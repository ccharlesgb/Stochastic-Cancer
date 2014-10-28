# -*- coding: utf-8 -*-
"""
Created on Thu Oct 23 15:22:10 2014

@author: Jonny
"""
import TwoSpecies

#Initialize the Gillespie simulator with 10 cells
N=100
mySim = TwoSpecies.Gillespie(N)
dataPointCount = 25
mySim.j = int(N/2)
mySim.r0=1.0
mySim.r1=1.0


def GetGammaj(sim):
    return sim.GetTJminus()/sim.GetTJplus()

def GetFixProb1(sim):
    sumf=0.0
    for k in range(1,N):
        prod=1.0    
        for j in range(1,k+1):
            sim.j=j
            prod*=GetGammaj(sim)
        sumf+=prod
    return 1/(1+sumf)

def GetFixTime1(sim):
    sum_k=0    #initialise sum over the upper limit of the the sum over TJplus 
    for k in range(1,N):    
        sum_Tplus=0 #initialise sum over TJplus
        for l in range(1,k+1):
            sim.j=l           #at each loop, vary j, so that gamma_j is different
            sum_Tplus+=(1/sim.GetTJplus())
            prod_of_gamma=1    #initialise the product over gamma j
            for j in range(l+1,k+1):
                sim.j=j        
                prod_of_gamma*=GetGammaj(sim)
            sum_Tplus*=prod_of_gamma
        sum_k+=sum_Tplus 
    return (1/(GetFixProb1(sim)))*sum_k
    
def GetFixTimeJ(sim,j):
    sumf=0.0
    #First Term
    for k in range(j,N):
        prod=1.0    
        for m in range(1,k+1):
            sim.j=m
            prod*=GetGammaj(sim)
        sumf+=prod
    
    part1=-GetFixTime1(sim)*sumf
    
    sum1=0    
    for k in range(j,N):    
        sum2=0
        for l in range(1,k+1):
            sim.j=l
            sum2+=(1/sim.GetTJplus())
            prod=1    
            for m in range(l+1,k+1):
                sim.j=m        
                prod*=GetGammaj(sim)
            sum2*=prod
        sum1+=sum2
    part2=sum1
    
    return part1+part2
    