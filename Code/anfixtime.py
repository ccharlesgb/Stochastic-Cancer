# -*- coding: utf-8 -*-
"""
Created on Thu Oct 23 15:22:10 2014

@author: Jonny
"""
import TwoSpecies
import matplotlib.pyplot as plt


def GetGammaj(sim):
    return sim.r1/sim.r0 #simplified version of gamma

def GetFixProb1(sim):
    sumf=0.0 #intialise sum
    for k in range(1,sim.N): #loop over 1
        prod=1.0    
        for j in range(1,k+1):
            prod*=GetGammaj(sim) 
        sumf+=prod

    #sim.j=1.0
    #print("FP1 for {0} is {1}".format(sim.j, 1.0/(1.0+sumf)))
    return 1.0/(1.0+sumf)

def GetFixTime1(sim):
    
    sum_k=0    #initialise sum over the upper limit of the the sum over TJplus 
    for k in range(1,sim.N):    
        sum_Tplus=0.0 #initialise sum over TJplus
        
        for l in range(1,k+1):
            #sim.j=l           #at each loop, vary j, so that gamma_j is different
            sum_Tplus_term = (1.0/sim.GetTJplus(l))
            #print("J={0} TJplus={1}".format(l,sim.GetTJplus()))
            
            
            prod_of_gamma=1.0    #initialise the product over gamma j            
            for j in range(l+1,k+1):
                #sim.j=j        
                prod_of_gamma*=GetGammaj(sim)
        
           
            sum_Tplus_term *= prod_of_gamma
            sum_Tplus += sum_Tplus_term
            
        sum_k+=sum_Tplus
    
    print(sum_k)
    #print("fixtime1 ={0}".format((1.0/(GetFixProb1(sim)))*sum_k))
    return (GetFixProb1(sim))*sum_k
    #return (1.0/0.1)*sum_k
    
def GetFixTimeJ(sim,j,fixTime1):
    sumf=0.0 #initialise sum over the upper limit of the the product of gamma j 
    #First Term
    for k in range(j, sim.N): #loop over the upper limit of the product of gamma j
        prod=1.0     #initialise product
        
        for m in range(1,k+1):  #loop over gamma j product
            #sim.j=m             #on each loop vary j so that gamma j is different
            prod*=GetGammaj(sim) #multiply by gamma j
            
        sumf+=prod        #multiply the product to each term in the sum
    
    #print("FixTime1 {0}".format(GetFixTime1(sim)))
    part1=-fixTime1*sumf #multiply by the fixation time of j=1 to get the first term
    #part1=-GetFixTime1(sim)*sumf #multiply by the fixation time of j=1 to get the first term
    #-GetFixTime1(sim)*sumf
    
    #second term
    
    sum1=0.0    #initialise first sum in second term
    
    for k in range(j,sim.N):    #loop upper index of subsequent sum and product
        sum2=0.0     #initialise second sum in second term
        
        for l in range(1,k+1): #loop over j's to get sum of 1/Tj+
            #sim.j=l         #change j every time
            sum2_term =(1.0/sim.GetTJplus(l))  #perform summation
            #print("J={0} TJplus={1}".format(l,sim.GetTJplus()))
            prod=1.0       #intialise product
            
            for m in range(l+1,k+1): #loop over product of gamma j for l+1 to k+1
                #sim.j=m        
                prod*=GetGammaj(sim)
                
            sum2_term *=prod #add the product of gamma j to each term in the series
            
            sum2 += sum2_term
        sum1+=sum2
        
    part2=sum1
    print("P1 {0} P2 {1}".format(part1,part2))
    return part1+part2
    
#Initialize the Gillespie simulator with 10 cells
mySim = TwoSpecies.Gillespie(100)
mySim.timeLimit = 1000000
dataPointCount = 100

mySim.r0=0.5
mySim.r1=0.5

dataPointsX = []
dataPointsY = []
fixTime1 = GetFixTime1(mySim)
for i in range(1,mySim.N):
    dataPointsX.append(i)
    dataPointsY.append(GetFixTimeJ(mySim,i,fixTime1))

plt.plot(dataPointsX, dataPointsY, label = "Theoretical")
plt.xlabel("Initial Number of Cell Types - j")
plt.ylabel("Theoretical Fixation Time - t_j")
plt.title("Theoretical Fixation Time for {0} Cells".format(mySim.N))