# -*- coding: utf-8 -*-
"""
Created on Thu Oct 23 15:22:10 2014

@author: Jonny
"""
#import TwoSpecies
import math

def GetGammaj(sim, param):
#    return sim.r1/sim.r0 #simplified version of gamma
    return param.GetTIJ(1,0, sim.n)/param.GetTIJ(0,1,sim.n)
    
def GetFixProb1(sim, param):
    sumf=0.0 #intialise sum
    for k in range(1,sim.N): #loop over 1
        prod=1.0    
        for j in range(1,k+1):
            prod*=GetGammaj(sim) 
        sumf+=prod

    #sim.j=1.0
    #print("FP1 for {0} is {1}".format(sim.j, 1.0/(1.0+sumf)))
    return 1.0/(1.0+sumf)
    
def GetFixTimeJ(sim, param,j):
    sumf=0.0 #initialise sum over the upper limit of the the product of gamma j 
    #First Term
    for k in range(j, param.N): #loop over the upper limit of the product of gamma j
        prod=1.0     #initialise product
        
        for m in range(1,k+1):  #loop over gamma j product
            sim.n[1] = m             #on each loop vary j so that gamma j is different
            prod*=GetGammaj(sim, param) #multiply by gamma j
            
        sumf+=prod        #multiply the product to each term in the sum

    Rate = param.u[0] * param.N
    
    part1 = 1.0/Rate * sumf
    
    #second term
    
    sum1=0.0    #initialise first sum in second term
    
    for k in range(j,param.N):    #loop upper index of subsequent sum and product
        sum2=0.0     #initialise second sum in second term
        
        for l in range(1,k+1): #loop over j's to get sum of 1/Tj+
            sim.n[1] = l         #change j every time
            param.PreSim(sim)
            sum2_term = 1.0/(param.GetTIJ(0,1, sim.n) / param.N)  #perform summation
            #print("J={0} TJplus={1}".format(l,sim.GetTJplus()))
            prod=1.0       #intialise product
            
            for m in range(l+1,k+1): #loop over product of gamma j for l+1 to k+1
                sim.n[1] = m        
                prod *= GetGammaj(sim, param)
                
            sum2_term *=prod #add the product of gamma j to each term in the series
            
            sum2 += sum2_term
        sum1+=sum2
        
    part2=sum1 / param.N
    #print("P1 {0} P2 {1}".format(part1,part2))
    return part1+part2