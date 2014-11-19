import random
import math

#Class for simulating cancer dynamics with the Gillespie algorithm
class Gillespie:
    def __init__(self, numCells):
        #Sim time parameters
        self.curTime = 0.0
        self.timeLimit = 100.0
        self.simSteps = 0
        
        #initialise population history
        self.populationHistory = 0      
        self.tHist = []
        self.jHist = []
        self.minusjHist = []
        
        #The fitness
        self.r0 = 1.0
        self.r1 = 2.0
        
        #the mutation rate
        self.u1=0.1

        self.N = numCells

        #Type Counts
        self.ij=0      
        self.j = self.ij
  

        self.lambd = self.GetLambda()
        
        self.preSim = 0
    
    #Reset parameters to default values
    def ResetSim(self):
        self.curTime = 0.0
        self.simSteps = 0
        self.j = self.ij
        
        self.tHist = []
        self.jHist = []
        self.minusjHist = []

    #Helper function to renormalize fitness of current cell population
    def AvgFitness(self):
        return (self.r1*self.j + self.r0*(self.N-self.j))
    
    #Reaction probability for j -> j+1
    def GetTJplus(self):
        top = (self.N-self.j)*( (self.j * self.r1) + self.u1 * self.r0 * (self.N - self.j))
        return top / (self.AvgFitness()*self.N)
    
    #probability for j -> j-1
    def GetTJminus(self):
        top = self.j * (self.N - self.j) * self.r0 * (1.0 - self.u1)
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
        return 0
    
    def RecordFrame(self):
        self.tHist.append(self.curTime)
        self.jHist.append(self.j)

    #Simulate loop
    #Uses standard gillespie algorithm and chooses event until fixated or out of time
    def Simulate(self):
        self.ResetSim()
        while self.curTime < self.timeLimit:
            if self.preSim:
                self.preSim(self)
            if self.populationHistory == 1:
                self.RecordFrame()
            if self.Fixated(): #Are we at absorbing state? If so quit
                return
            timestep = self.GetTimeStep() #How much time until the next event?
            self.ChooseEvent() #Chose what kind of event and update cell counts

            self.curTime += timestep #Increment time
            self.simSteps+= 1 #Increase event count
#END CLASS DEFINTION


#Pulse wave function
def PulseWave(time, amp, frequency, width_frac=0.5, offset=0):
    if frequency == 0:
        return amp #no frequency means that there is no treatment - therefore fitness/mutation rate should not be varied from its maxium
    
    period=1.0/frequency
    time+=offset*period
    width=period*width_frac
    x,y = divmod(time, float(period))
    if (y < width):
        return amp
    else:
        return 0

'''
#Pulse wave function
def PulseWave(time, amp, width, period):
    x,y = divmod(time, float(period))
    if (y < width):
        return amp
    else:
        return 0.0
'''