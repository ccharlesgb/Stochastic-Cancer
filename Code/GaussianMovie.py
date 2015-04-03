"""
A simple example of an animated plot
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import wright_fisher
import TauLeapParam
import matplotlib.pyplot as plt
import math

def GetTauJ(param):
    summation = 0.0
    s = param.r[1] - param.r[0]
    for i in range(0, param.cellTypes):
        top = math.log(s / (param.u[i] * param.d))
        bottom = 2.0 * s * math.log(param.N)
        summation += (top / bottom)
    return summation;


cellTypes = 21
population = 1e9

myWF = wright_fisher.wright_fisher(cellTypes)
myHist = TauLeapParam.Hist(cellTypes)
myParam = TauLeapParam.Params(cellTypes)

myWF.stopAtAppear = 1
myWF.history = myHist
myWF.stepLimit = 10000
myWF.useApproxTheta = 0
myWF.params = myParam

myParam.n0[0] = population
myParam.u = [1e-7] * cellTypes
myParam.d = 100
myParam.uNotConst = 0

s = 1e-2
d = 1.8

start = 0.1
end = 1e-2

sum_s_j = 0.0
for i in range(0,cellTypes):
    myParam.r[i] = 1.0 + s * (float(i) / math.sqrt(cellTypes))*(float(i) / math.sqrt(cellTypes))
    s_j = 0.1 - ((start - end)/((cellTypes-1)**2))*(i**2)
    sum_s_j += s_j
    myParam.r[i] = 1.0 + sum_s_j
    myParam.r[i] = math.pow(1.0 + s, i)

myWF.Simulate()

def GetDistribution(t):
    dataY = []
    for i in range(0, cellTypes): 
        dataY.append(myHist.histArray[i][t])
    return dataY

fig, ax = plt.subplots(nrows = 2)

#x = np.arange(0.0, cellTypes, 1.0)
x = []
sj_x = []
for i in range(0, cellTypes):
    sj_x.append(myParam.r[i] - 1.0)
    x.append(i)
print(x)
print(GetDistribution(0))

line, = ax[0].semilogy(x, GetDistribution(0), 'o-')
ax[0].set_ylim(10^9, 0)
line2, = ax[1].semilogy(sj_x, GetDistribution(0), 'o-')
ax[1].set_ylim(10^9, 0)

ax[0].set_xlabel("j")
ax[1].set_xlabel("s_j")

ax[0].set_ylabel("N_j")
ax[1].set_ylabel("N_sj")

#x = np.arange(0, 2*np.pi, 0.01)        # x-array
#line, = ax.plot(x, np.sin(x))

def animate(i):
    yDat = GetDistribution(i)
    line.set_ydata(yDat)  # update the data
    line2.set_ydata(yDat)
    return line,line2,

#Init only required for blitting to give a clean slate.
def init():
    line.set_ydata(np.ma.array(x, mask=True))
    line2.set_ydata(np.ma.array(sj_x, mask=True))
    return line,line2,


ani = animation.FuncAnimation(fig, animate, np.arange(0, myWF.curTime,10), init_func=init,
    interval=50, blit=True)
plt.show()