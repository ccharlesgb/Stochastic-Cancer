# -*- coding: utf-8 -*-
"""
Created on Thu Oct 23 15:22:10 2014

@author: Jonny
"""
import TwoSpecies

#Initialize the Gillespie simulator with 10 cells
N=10
mySim = TwoSpecies.Gillespie(N)
dataPointCount = 25
mySim.j = int(N/2)
mySim.r0=1.0
mySim.r1=2.0


def GetGammaj(sim):
    return sim.r1/sim.r0 #simplified version of gamma

def GetFixProb1(sim):
    sumf=0.0 #intialise sum
    for k in range(1,N-1): #loop over 1
        prod=1.0    
        for j in range(1,k+1):
            sim.j=j #change j 
            prod*=GetGammaj(sim) 
        sumf+=prod
    return 1/(1+sumf)

def GetFixTime1(sim):
    sum_k=0    #initialise sum over the upper limit of the the sum over TJplus 
    for k in range(1,N-1):    
        sum_Tplus=0 #initialise sum over TJplus
        for l in range(1,k+1):
            sim.j=l           #at each loop, vary j, so that gamma_j is different
            print(l)            
            sum_Tplus+=(1/sim.GetTJplus())
            prod_of_gamma=1    #initialise the product over gamma j
            for j in range(l+1,k+1):
                sim.j=j        
                prod_of_gamma*=GetGammaj(sim)
            sum_Tplus*=prod_of_gamma
        sum_k+=sum_Tplus 
    return (1/(GetFixProb1(sim)))*sum_k
    
def GetFixTimeJ(sim,j):
    sumf=0.0 #initialise sum over the upper limit of the the product of gamma j 
    #First Term
    for k in range(j,N-1): #loop over the upper limit of the product of gamma j
        prod=1.0     #initialise product
        for m in range(1,k+1):  #loop over gamma j product
            sim.j=m             #on each loop vary j so that gamma j is different
            prod*=GetGammaj(sim) #multiply by gamma j
            print(sim.j)
            print(GetGammaj(sim))
        sumf*=prod        #multiply the product to each term in the sum
    
    part1=-GetFixTime1(sim)*sumf #multiply by the fixation time of j=1 to get the first term
    
    sum1=0    #initialise first sum in second term
    for k in range(j,N-1):    #loop upper index of subsequent sum and product
        sum2=0                  #initialise second sum in second term
        for l in range(1,k+1): #loop over j's to get sum of 1/Tj+
            sim.j=l         #change j every time
            sum2+=(1/sim.GetTJplus())  #perform summation
            prod=1       #intialise product
            for m in range(l+1,k+1): #loop over product of gamma j for l+1 to k+1
                sim.j=m        
                prod*=GetGammaj(sim)
            sum2*=prod #add the product of gamma j to each term in the series
        sum1+=sum2
    part2=sum1
    
    return part1+part2
    