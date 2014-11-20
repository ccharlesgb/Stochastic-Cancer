# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 10:52:59 2014

@author: Jonny
"""
import TwoSpecies
import anfixtime

N=10 #define number of cells
mySim = TwoSpecies.Gillespie(N) #initialise sim
mySim.r0 = 1.0 #set system parameters (for now)
mySim.r1 = 1.0
mySim.u1 = 0.01

r1_max_pulse=1.3
r1_min_pulse=0.7

#set up helper function to find w_k
def wk(sim,k):
    #work out the first term
    prod=1.0
    for i in range(1,k):
        sim.j=i        
      
        term=anfixtime.GetGammaj(sim)
        prod*=term
    term1=prod/(sim.u1)
    
    #work out the second term
    sum_k=0.0    
    for j in range(1,k):
        sim.j=j        
        prod_l = 1.0
        for l in range(j+1,k):
            term=anfixtime.GetGammaj(sim)
            prod_l*=term
        sum_k+=(1.0/sim.GetTJplus() )*prod_l
    term2=-sum_k
    
    return term1 + term2
def theorypulse(sim, max_r1,min_r1, frequency): 
    global N
    sim.j=0   
    k=1
    T=1.0/(2.0*frequency)    
    #T=10.0
    timecount=0.0
    treatment=0
    while k<N:
        if treatment == 0:
            sim.r1=r1_max_pulse
        if treatment == 1:
            sim.r1=r1_min_pulse
        sum_wk=0.0
        while sum_wk > -T and k < N:
            sum_wk+=wk(sim,k)
            k+=1
        print(k)
        timecount+=T
        print("treatment:{0}".format(treatment))        
        treatment = not treatment
    return timecount
    

quicktest=theorypulse(mySim, 1.3,0.7,0.025)
print("Fix time for frequency 0.025 = {0}".format(quicktest))

