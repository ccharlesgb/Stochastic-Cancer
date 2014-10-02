import random
import math
import matplotlib.pyplot as plt

random.seed()

#Struct for passing simulation parameters
class SimParam:
    def __init__(self, numCells):
        self.curTime = 0.0
        self.timeLimit = 100.0
        
        #The fitness
        self.r0 = 1.0
        self.r1 = 2.0
        self.r2 = 1.0

        #Mutation Rates
        self.u1 = 0.1
        self.u2 = 0.1

        self.N = numCells

        #Type Counts
        self.n0 = self.N
        self.n1 = 0
        self.n2 = 0

        self.cells = [] #Is this actually needed?
        for i in range(0, self.N):
            self.cells.insert(i,0)

        self.lambd = self.GetLambda()

    def ResetSim(self):
        self.curTime = 0.0
        self.n0 = self.N
        self.n1 = 0
        self.n2 = 0
        self.cells = []
        for i in range(0, self.N):
            self.cells.insert(i,0)
            
    def AvgFitness(self):
        return (self.r0*self.n0 + self.r1*self.n1 + self.r2 * self.n2)/self.N
        
    #Reaction probability for cell from 1->0
    def GetT10(self):
        top = (1 - self.u1) * self.r0 * (float(self.n0)/self.N) * self.n1
        return top / self.AvgFitness()

    def GetT20(self):
        top = (1 - self.u1) * self.r0 * (float(self.n0)/self.N) * self.n2
        return top / self.AvgFitness()

    def GetT01(self):
        top = (self.u1 * self.r0 * (float(self.n0)/self.N) + (1-self.u2) * self.r1 * (float(self.n1)/self.N)) * self.n0
        return top / self.AvgFitness()

    def GetT21(self):
        top = (self.u1 * self.r0 * (float(self.n0)/self.N) + (1-self.u2) * self.r1 * (float(self.n1)/self.N)) * self.n2
        return top / self.AvgFitness()

    def GetT02(self):
        top = (self.u2 * self.r1 * (float(self.n1)/self.N) + self.r2 * (float(self.n2)/self.N)) * self.n0
        return top / self.AvgFitness()

    def GetT12(self):
        top = (self.u2 * self.r1 * (float(self.n1)/self.N) + self.r2 * (float(self.n2)/self.N)) * self.n1
        return top / self.AvgFitness()

    #Exponential parameter for frequency of events
    def GetLambda(self):
        self.lambd = self.GetT10() + self.GetT20() + self.GetT01() + self.GetT21() + self.GetT02() + self.GetT12()
        return self.lambd

    #Returns an exponentially distributed number based on the lambda parameter
    def GetTimeStep(self):
        return 1.0/(self.GetLambda() * math.log(1.0/random.random()))

    #Chose and execute which event to carry out. Updates cell array and population counts
    def ChooseEvent(self):
        rand = random.random()
        lam = self.GetLambda()
        #Cell Type 1->0
        threshold = self.GetT10()/lam
        #print(threshold)
        if (rand < threshold):
            index = self.GetFirstIndexOf(1)
            self.cells[index] = 0
            self.n0 += 1
            self.n1 -= 1
            return
        #Cell Type 2->0
        threshold += self.GetT20()/lam
        #print(threshold)
        if (rand < threshold):
            index = self.GetFirstIndexOf(2)
            self.cells[index] = 0
            self.n0 += 1
            self.n2 -= 1
            return
        #Cell Type 0->1
        threshold += self.GetT01()/lam
        #print(threshold)
        if (rand < threshold): 
            index = self.GetFirstIndexOf(0)
            self.cells[index] = 1
            self.n1 += 1
            self.n0 -= 1
            return
        #Cell Type 2->1
        threshold += self.GetT21()/lam
        #print(threshold)
        if (rand < threshold):
            index = self.GetFirstIndexOf(2)
            self.cells[index] = 1
            self.n1 += 1
            self.n2 -= 1
            return
        #Cell Type 0->2
        threshold += self.GetT02()/lam
        #print(threshold)
        if (rand < threshold): 
            index = self.GetFirstIndexOf(0)
            self.cells[index] = 2
            self.n2 += 1
            self.n0 -= 1
            return
        #Cell Type 1->2
        threshold += self.GetT12()/lam
        #print(threshold)
        if (rand < threshold): 
            index = self.GetFirstIndexOf(1)
            self.cells[index] = 2
            self.n2 += 1
            self.n1 -= 1
            return

    #Helper function finds the first instance of a type of cell in the cell array
    def GetFirstIndexOf(self, celltype):
        for i in range(0, self.N):
            if self.cells[i] == celltype:
                return i

    def Fixated(self):
        return self.n2 >= self.N
    
    def Simulate(self):
        self.ResetSim()
        #print("n0 = %i" % mySim.n0)
        #print("n1 = %i" % mySim.n1)
        #print("n2 = %i" % mySim.n2)
        while self.curTime < self.timeLimit:
            timestep = self.GetTimeStep()
            self.ChooseEvent()
            #print("Event:")
            #print("n0 = %i" % mySim.n0)
            #print("n1 = %i" % mySim.n1)
            #print("n2 = %i" % mySim.n2)
            if self.Fixated():
                return
            self.curTime += timestep
#END CLASS DEFINTION


mySim = SimParam(10)
mySim.timeLimit = 1000
mySim.u2 = 0.001
dataPointCount = 25

simsPerDataPoint = 100

dataPointsX = []
dataPointsY = []
for i in range(0, dataPointCount):
    dataPointsX.insert(i,0.0)
    dataPointsY.insert(i,0.0)

for curPoint in range(0, dataPointCount):
    fixationCounts = 0
    mySim.r1 = 0.2 + (3.0-0.2) * float(curPoint)/float(dataPointCount) #sweep type-1 fitness from 0.2 to 3.0
    print("Type-1 Fitness: %f" % mySim.r1)
    for sim in range(0, simsPerDataPoint):
        mySim.Simulate()

        #print("Simulation Complete")
        #print("n0 = %i" % n0)
        #print("n1 = %i" % n1)
        #print("n2 = %i" % n2)

        if mySim.n2 == mySim.N:
            fixationCounts += 1
    dataPointsX[curPoint] = mySim.r1
    dataPointsY[curPoint] = float(fixationCounts) / float(simsPerDataPoint)

for i in range(0, dataPointCount):
    print("r1: {0} Fixation: {1}".format(dataPointsX[i],dataPointsY[i]))
    
filename = "sim_N={0}_r0={1}_r1={2}_r2={3}_u1={4}_u2={5}_SPDP={6}".format(mySim.N, mySim.r0, mySim.r1, mySim.r2, mySim.u1, mySim.u2, simsPerDataPoint)

f = open(filename, 'w')

for i in range(0, dataPointCount):
    f.write(str(dataPointsX[i]) + "    " + str(dataPointsY[i]) + "\n")

f.close()

plt.plot(dataPointsX, dataPointsY)
plt.xlabel("Type 1 Fitness r1")
plt.ylabel("Type 2 Fixation %")
plt.show()




        
