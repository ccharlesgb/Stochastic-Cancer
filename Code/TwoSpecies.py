import random
import math

random.seed()

class SimParam:
    def __init__(self, numCells):
        self.timeLimit = 0.0
        
        self.r0 = 0.0
        self.r1 = 0.0
        
        self.u1 = 0.0
        
        self.N = 0
        
        self.j = 0

        
        

#Class for simulating cancer dynamics with the Gillespie algorithm
class Gillespie:
    def __init__(self, numCells):
        #Sim time parameters
        self.curTime = 0.0
        self.timeLimit = 100.0
        self.simSteps = 0

        
        #The fitness
        self.r0 = 1.0
        self.r1 = 2.0


        self.N = numCells

        #Type Counts
        self.j = self.N
  

        self.lambd = self.GetLambda()
        
        self.preSim = 0
    
    #Reset parameters to default values
    def ResetSim(self):
        self.curTime = 0.0
        self.simSteps = 0
        self.j = (self.N)/2
        
    
    #Helper function to renormalize fitness of current cell population
    def AvgFitness(self):
        return (self.r0*self.j + self.r1*(self.N-self.j))
        
    #Reaction probability for cell from 1->0
    def GetTJplus(self):
        top = self.j*self.r0 *(self.N-self.j)
        print("TJplus = {0}".format(top / (self.AvgFitness()*self.N)))
        return top / (self.AvgFitness()*self.N)

    def GetTJminus(self):
        top = (self.N-self.j)*self.r1*self.j
        print("TJminus {0}".format(top / (self.AvgFitness()*self.N)) )        
        return top / (self.AvgFitness()*self.N)


    #Exponential parameter for frequency of events
    def GetLambda(self):
        self.lambd = self.GetTJplus() + self.GetTJminus()    
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
        #Cell Type j->j+1
        threshold = self.GetTJplus()/lam
        #print(threshold)
        if (rand < threshold):
            self.j += 1
            return
        #Cell Type j->j-1
        threshold += self.GetTJminus()/lam
        #print(threshold)
        if (rand < threshold):
            self.j -= 1
            return
        
    
    #Returns true if we are at the absorbing state of j == N
    def Fixated(self):
        if self.j>=self.N:
            return 1
        if self.j<=0:
            return 1
        return 0
        

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

#Pulse wave function
def PulseWave(time, amp, width, wavelength):
    x,y = divmod(time, wavelength)
    if (y < width):
        return amp
    else:
        return 0.0

