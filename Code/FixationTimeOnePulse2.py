# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import anfixtime
import TwoSpecies
import time
import math

#THEORY

def Getwk(sim, k):
    #First Term
    prod1 = 1.0
    for i in range(1, k):
        sim.j = i
        prod1 = prod1 * anfixtime.GetGammaj(sim)
        
    part1 = (1.0/(sim.u1 * sim.N)) * prod1
    #print("WK INITIAL {0}".format(wk))
    #Second Term
    sum1 = 0.0
    for j in range(1, k):
        prodgamm = 1.0
        for l in range(j+1, k):
            sim.j = l
            prodgamm = prodgamm * anfixtime.GetGammaj(sim)
            
        sim.j = j
        sum1 = sum1 + (1.0/(sim.GetTJplus() / sim.N)) * prodgamm
    
    part2 = sum1 / sim.N
    
    #print("PART 1 {0} PART 2 {1}".format(part1, part2))
    
    return part1 + part2

def FindMeanJ(sim, t_mut):
    #print("T_MUT IS {0}".format(t_mut))
    mean_j = -1

    last_sum = -1.0    
    
    for ul in range(1, sim.N):
        sum_wk = 0.0
        for k in range(1, ul + 1):
            sum_wk += Getwk(sim, k)
        #print("SUM WK = {0}: {1}".format(ul, sum_wk))
        
        if sum_wk >= t_mut:
            if last_sum != -1.0 and (abs(last_sum - t_mut) < abs(sum_wk - t_mut)):
                print("STUFF")
                print(abs(last_sum - t_mut))
                print(abs(sum_wk - t_mut))
                mean_j = ul - 1
            else:
                mean_j = ul
            break
        last_sum = sum_wk
    if mean_j == -1: #Didnt find it
        mean_j = sim.N
    print("FOUND IT AT J = {0}".format(ul))
    return mean_j
    
def FixTimePulse(sim, t_mut, quit_u1, smoke_u1):
    #Get the mean j that we will be at in the fix time
    sim.u1 = smoke_u1
    mean_j = FindMeanJ(sim, t_mut)
    
    sim.u1  = quit_u1
    fix_from_mean = anfixtime.GetFixTimeJ(sim, mean_j)

    fixtime = t_mut + fix_from_mean
    
    return fixtime

#Initialize the Gillespie simulator with 10 cells
mySim = TwoSpecies.Gillespie(10)

mySim.ij = 0
mySim.N = 100
mySim.timeLimit = 100.0
mySim.r0 = 1.0
mySim.r1 = 1.0

minSmokeTime = 0.0
maxSmokeTime = 250

smokeArea_u1 = 1.0

quit_u1 = 0.01

smoke_u1 = quit_u1 * 10.0

min_mut_factor = 2.0
max_mut_factor = 20.0


#Sweep the parameter r1 from 0.2 to 3.0 and run many simulations per data point
#Gets an idea on how likely cancer fixation is to occur for this parameter
dataPointCount = 5

fixError = []

#Initialize the array with default values
dataPointsX = []
dataPointsY = []

dataPointsY2 = []

dataST = []
dataMU = []

for curPoint in range(0, dataPointCount):
    startTime = time.clock() #Algorithm benchmarking
    
    fixationCounts = 0
    totFixTime = 0.0
    print("Current Data Point = {0}/{1} ({2}%)".format(curPoint + 1, dataPointCount, 100.0 * float(curPoint+1.0)/dataPointCount))
    
    #mut_factor = (max_mut_factor - min_mut_factor) * ((dataPointCount) / float(curPoint + 1.0)) + min_mut_factor
    mut_factor = (max_mut_factor - min_mut_factor) * ((1.0) / float(curPoint + 1.0)) + min_mut_factor
    smoke_u1 = quit_u1 * mut_factor
    smokeTime = smokeArea_u1 / (smoke_u1 - quit_u1)
    
    dataST.append(smokeTime)
    dataMU.append(smoke_u1)
    
    print("MUT FACTOR: {0}".format(smoke_u1))
    print("T SMOKE: {0}".format(smokeTime))
    
    #print("CHECK DOSE: {0}".format((smoke_u1 - quit_u1) * smokeTime))
    
    dataPointsX.append(smokeTime)
    dataPointsY.append(FixTimePulse(mySim, smokeTime, quit_u1, smoke_u1))
    
    mySim.u1 = quit_u1
    dataPointsY2.append(anfixtime.GetFixTimeJ(mySim, 0))


plt.plot(dataPointsX, dataPointsY)
#plt.plot(dataPointsX, dataPointsY2, ':')
plt.xlabel("Mutagen Exposure Time")
plt.ylabel("Fixation Time")

'''

#Initialize the array with default values
dataPointsX = []
dataPointsY = []

dataPointsY2 = []

dataST = []
dataMU = []

#Pre Simulate callback (called every frame before a timestep)
def increaseMutation(sim):
    global curPoint
    global dataPointCount
    global smokeTime
  
    if sim.curTime < smokeTime: #Still smoking
        sim.u1 = smoke_u1
    else:
        sim.u1 = quit_u1

mySim.preSim = increaseMutation #IMPORTANT assign the callback (called in the class sim loop)

total_dist = []
raw_avgs = []

for curPoint in range(0, dataPointCount):
    startTime = time.clock() #Algorithm benchmarking
    
    fixationCounts = 0
    totFixTime = 0.0
    print("Current Data Point = {0}/{1} ({2}%)".format(curPoint + 1, dataPointCount, 100.0 * float(curPoint+1.0)/dataPointCount))
    
    #mut_factor = (max_mut_factor - min_mut_factor) * ((dataPointCount) / float(curPoint + 1.0)) + min_mut_factor
    mut_factor = (max_mut_factor - min_mut_factor) * ((1.0) / float(curPoint + 1.0)) + min_mut_factor
    smoke_u1 = quit_u1 * mut_factor
    smokeTime = smokeArea_u1 / (smoke_u1 - quit_u1)
    
    dataST.append(smokeTime)
    dataMU.append(smoke_u1)
    
    print("MUT FACTOR: {0}".format(smoke_u1))
    print("T SMOKE: {0}".format(smokeTime))
    
    print("CHECK DOSE: {0}".format((smoke_u1 - quit_u1) * smokeTime))
    
    
    dataPointsX.append(smokeTime)
    mySim.u1 = smoke_u1
    dataPointsY.append(FindMeanJ(mySim, smokeTime))
    
    indiv_dist = [] 

    sum_ind = 0.0
    
    for sid in range(0, 10000):
        mySim.timeLimit = smokeTime
        mySim.Simulate()
        indiv_dist.append(mySim.j)
        sum_ind += mySim.j
        
    raw_avgs.append(sum_ind / 10000)
    
    total_dist.append([])
    total_dist[curPoint] = indiv_dist


for i in range(0, dataPointCount):
    plt.figure()
    plt.hist(total_dist[i], 20)
    plt.title("MEAN J IS {0} (Exact: {1})".format(dataPointsY[i], raw_avgs[i]))
    plt.xlim(0, 100)
'''