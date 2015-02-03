import random
import math
import matplotlib.pyplot as plt

random.seed()

nx = 0.0
ny = 0.0

Omega = 100
b = 100

#Rate 1
def GetW1():
    return Omega;
    
#Rate 2
def GetW2():
    return nx;
    
#Rate 3
def GetW3():
    return b * nx
    
#Rate 4
def GetW4():
    return (nx*nx*ny)/(Omega*Omega);

#Exponential parameter for frequency of events
def GetLambda():
    lambd = GetW1() + GetW2() + GetW3() + GetW4()
    return lambd

#Returns an exponentially distributed number based on the lambda parameter
def GetTimeStep():
    return 1.0/GetLambda() * math.log(1.0/random.random())
    #return random.expovariate(self.GetLambda())

timeLimit = 30.0

curT = 0.0

tHist = []
nxHist = []
nyHist = []

def ChooseEvent():
    lambd = GetLambda()
    ran = random.random()
    global nx
    global ny
    threshold = GetW1()/lambd
    if (ran < threshold):
        nx += 1
        return
    threshold += GetW2()/lambd
    if (ran < threshold):
        nx -= 1
        return
    threshold += GetW3()/lambd
    if (ran < threshold):
        nx -= 1
        nx += 1
        return
    threshold += GetW4()/lambd
    if (ran < threshold):
        nx += 1
        ny -= 1 
        return

while curT < timeLimit:
    deltaT = GetTimeStep()
    
    ChooseEvent()

    curT += GetTimeStep()

    tHist.append(curT)
    nxHist.append(nx)
    nyHist.append(ny * 2.0)
    
plt.plot(tHist, nxHist)
plt.plot(tHist, nyHist)
plt.show()


