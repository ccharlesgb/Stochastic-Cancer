import random
import math

random.seed()

#Class for simulating cancer dynamics with the Gillespie algorithm
class Gillespie:
    def __init__(self, numCells):
        #Sim time parameters
        self.curTime = 0.0
        self.timeLimit = 100.0
        self.simSteps = 0
        
        self.populationHistory = 0        
        
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

        self.lambd = self.GetLambda()
        
        self.preSim = 0

    #Reset parameters to default values
    def ResetSim(self):
        self.curTime = 0.0
        self.simSteps = 0
        self.n0 = self.N
        self.n1 = 0
        self.n2 = 0
    
    #Helper function to renormalize fitness of current cell population
    def AvgFitness(self):
        return (self.r0*self.n0 + self.r1*self.n1 + self.r2 * self.n2)
        
    #Reaction probability for cell from 1->0
    def GetT10(self):
        top = (1 - self.u1) * self.r0 * float(self.n0) * self.n1
        return top / self.AvgFitness()

    def GetT20(self):
        top = (1 - self.u1) * self.r0 * float(self.n0) * self.n2
        return top / self.AvgFitness()

    def GetT01(self):
        top = (self.u1 * self.r0 * self.n0 + (1-self.u2) * self.r1 * self.n1) * self.n0
        return top / self.AvgFitness()

    def GetT21(self):
        top = (self.u1 * self.r0 * self.n0 + (1-self.u2) * self.r1 * self.n1) * self.n2
        return top / self.AvgFitness()

    def GetT02(self):
        top = (self.u2 * self.r1 * self.n1 + self.r2 * self.n2) * self.n0
        return top / self.AvgFitness()

    def GetT12(self):
        top = (self.u2 * self.r1 * self.n1 + self.r2 * self.n2) * self.n1
        return top / self.AvgFitness()

    #Exponential parameter for frequency of events
    def GetLambda(self):
        self.lambd = self.GetT10() + self.GetT20() + self.GetT01() + self.GetT21() + self.GetT02() + self.GetT12()
        return self.lambd

    #Returns an exponentially distributed number based on the lambda parameter
    def GetTimeStep(self):
        return 1.0/self.GetLambda() * math.log(1.0/random.random())
        #return random.expovariate(self.GetLambda())
        
    #Chose and execute which event to carry out. Updates population counts
    #Uses weighted random number between 0 and 1
    def ChooseEvent(self):
        rand = random.random()
        lam = self.GetLambda()
        #Cell Type 1->0
        threshold = self.GetT10()/lam
        #print(threshold)
        if (rand < threshold):
            self.n0 += 1
            self.n1 -= 1
            return
        #Cell Type 2->0
        threshold += self.GetT20()/lam
        #print(threshold)
        if (rand < threshold):
            self.n0 += 1
            self.n2 -= 1
            return
        #Cell Type 0->1
        threshold += self.GetT01()/lam
        #print(threshold)
        if (rand < threshold): 
            self.n1 += 1
            self.n0 -= 1
            return
        #Cell Type 2->1
        threshold += self.GetT21()/lam
        #print(threshold)
        if (rand < threshold):
            self.n1 += 1
            self.n2 -= 1
            return
        #Cell Type 0->2
        threshold += self.GetT02()/lam
        #print(threshold)
        if (rand < threshold): 
            self.n2 += 1
            self.n0 -= 1
            return
        #Cell Type 1->2
        threshold += self.GetT12()/lam
        #print(threshold)
        if (rand < threshold): 
            self.n2 += 1
            self.n1 -= 1
            return
    
    #Returns true if we are at the absorbing state of n2 == N
    def Fixated(self):
        return self.n2 >= self.N
        
    #Simulate loop
    #Uses standard gillespie algorithm and chooses event until fixated or out of time
    def Simulate(self):
        self.ResetSim()
        while self.curTime < self.timeLimit:
            if self.preSim:
                self.preSim(self)
            timestep = self.GetTimeStep() #How much time until the next event?
            self.ChooseEvent() #Chose what kind of event and update cell counts
            if self.Fixated(): #Are we at absorbing state? If so quit
                return
            self.curTime += timestep #Increment time
            self.simSteps+= 1 #Increase event count
#END CLASS DEFINTION

#Helper functions

#Saves an array of xy data to a file
def SaveXYToFile(filename, dataX, dataY):
    f = open(filename, 'w')

    #Dump data to file
    for i in range(0, len(dataX)):
        f.write(str(dataX[i]) + "    " + str(dataY[i]) + "\n")
    
    f.close()


