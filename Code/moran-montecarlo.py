import random

random.seed()

#Returns a cell type to divide based on the current population
def CellDivide():
    #Calculate Division Probabilities for each type
    totalFitness = r0*n0 + r1*n1 + r2*n2
    divProb0 = r0*n0 * (1-u1)
    divProb1 = r0*n0*u1 + r1*n1*(1-u2)
    divProb2 = r1*n1*u2 + r2*n2
    divProb0 /= totalFitness
    divProb1 /= totalFitness
    divProb2 /= totalFitness

    #print(divProb0)
    #print(divProb1)
    #print(divProb0 + divProb1)

    unitRand = random.random()
    if unitRand < divProb0:
        return 0
    if unitRand < divProb0 + divProb1:
        return 1
    if unitRand < divProb0 + divProb1 + divProb2:
        return 2
#

#Increment/Decrement the populations based on what cell died and which daughte cell was produced
def UpdatePopulations(new, old):
    if new == old: #Get out clause if the same type died as divided then nothing happened
        return

    global n0
    global n1
    global n2
    
    if new == 0:
        n0 += 1
    if new == 1:
        n1 += 1
    if new == 2:
        n2 += 1

    if old == 0:
        n0 -= 1
    if old == 1:
        n1 -= 1
    if old == 2:
        n2 -= 1
#

#Tick over timesteps and divide and kill a cell in each step
def Simulate():
    #MAIN SIM LOOP
    global N
    global n0
    global n1
    global n2
    #Initiate the Simulation
    cells = []
    for i in range(0, N):
        cells.insert(i,startType)
    n0 = N
    n1 = 0
    n2 = 0
    
    #Run the timesteps
    for t in range(0,stepCount):
        newCell = CellDivide()
        deathIndex = random.randint(0, N-1)
        UpdatePopulations(newCell, cells[deathIndex])
        cells[deathIndex] = newCell
#

#The fitness
r0 = 1.0
r1 = 2.0
r2 = 2.0

#Mutation Rates
u1 = 0.1
u2 = 0.001

#Type Counts
n0 = 0
n1 = 0
n2 = 0

cells = [];

startType = 0

#N = input("Enter number of cells:")
N = 10

for i in range(0, N):
    cells.insert(i,startType)

#stepCount = input("How many steps:")
stepCount = 100 * N #Is this correct?

dataPointCount = 25

simsPerDataPoint = 1000

dataPointsX = []
dataPointsY = []
for i in range(0, dataPointCount):
    dataPointsX.insert(i,0.0)
    dataPointsY.insert(i,0.0)

for curPoint in range(0, dataPointCount):
    fixationCounts = 0
    r1 = 0.2 + (3.0-0.2) * float(curPoint)/float(dataPointCount) #sweep type-1 fitness from 0.2 to 3.0
    print("Type-1 Fitness: %f" % r1)
    for sim in range(0, simsPerDataPoint):
        Simulate()

        #print("Simulation Complete")
        #print("n0 = %i" % n0)
        #print("n1 = %i" % n1)
        #print("n2 = %i" % n2)

        if n2 == N:
            fixationCounts += 1
    dataPointsX[curPoint] = r1
    dataPointsY[curPoint] = float(fixationCounts) / float(simsPerDataPoint)

for i in range(0, dataPointCount):
    print("r1: {0} Fixation: {1}".format(dataPointsX[i],dataPointsY[i]))
    
f = open('sim_output.dat', 'w')

for i in range(0, dataPointCount):
    f.write(str(dataPointsX[i]) + "    " + str(dataPointsY[i]) + "\n")

f.close()




        
