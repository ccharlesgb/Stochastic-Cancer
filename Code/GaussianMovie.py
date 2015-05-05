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

cellTypes = 21
population = 1e9

myWF = wright_fisher.wright_fisher(cellTypes)
myHist = TauLeapParam.Hist(cellTypes)
myParam = TauLeapParam.Params(cellTypes)

myWF.stopAtAppear = 1
myWF.history = myHist
myWF.stepLimit = 20000
myWF.useApproxTheta = 0
myWF.params = myParam

myParam.n0[0] = population
myParam.d = 100
myParam.uNotConst = 0

#Sinusoidal mutation rate
f = 0.25
for i in range(1, cellTypes):
    u_log = -7.0 + 3.0*math.cos(i*f * 2.0 * math.pi)
    u = math.pow(10.0, u_log)
    myParam.SetU(i,u)

'''
#Exponential mutation rate
for i in range(1, cellTypes):
    u_log = min(-8.0 + math.exp(i/20.0), -4.0)
    u = math.pow(10.0, u_log)
    myParam.SetU(i,u)
'''   
'''
#Switch
u_before = 1e-6
u_after = 1e-9
switch = 10
for i in range(1, cellTypes):
    if (switch + 1) > i:
        myParam.SetU(i, u_before)
    else:
        myParam.SetU(i, u_after)
'''
s = 1e-2
d = 1.8

start = 0.1
end = 1e-2

sum_s_j = 0.0
for i in range(0,cellTypes):
    myParam.r[i] = 1.0 + s * (float(i) / math.sqrt(cellTypes))*(float(i) / math.sqrt(cellTypes))
    #s_j = 0.1 - ((start - end)/((cellTypes-1)**2))*(i**2)
    #sum_s_j += s_j
    #myParam.r[i] = 1.0 + sum_s_j
    myParam.r[i] = math.pow(1.0 + s, i)

myParam.SetCompoundFitness(1e-2)
myParam.SetUAll(1e-7)

myWF.Simulate()

def GetDistribution(t):
    dataY = []
    for i in range(0, cellTypes): 
        dataY.append(myHist.histArray[i][t])
    return dataY
    
def GetAvgJ(t):
    avg = 0.0
    for i in range(0, cellTypes):
        avg = avg + (i * float(myHist.histArray[i][t])/myParam.N)
    return avg

def GetVelocityJ(t):
    if t == myWF.curTime:
        t = t - 1
    avg_t = GetAvgJ(t)
    avg_t1 = GetAvgJ(t+1)
    return avg_t1 - avg_t
    

fig, ax = plt.subplots(figsize = [16,9],facecolor='white')
#x = np.arange(0.0, cellTypes, 1.0)
x = []
#sj_x = []
for i in range(0, cellTypes):
    #sj_x.append(myParam.r[i] - 1.0)
    x.append(i)
print(x)
print(GetDistribution(0))

line, = ax.semilogy(x, GetDistribution(0), 'o-', label = "x_j(t) Distribution")
line2, = ax.semilogy([GetAvgJ(0)]*2, [1, 1e9], label = "Average j")
line3, = ax.semilogy([GetVelocityJ(0)]*2, [1, 1e9], label = "Velocity of avg j")
ax.set_ylim(0, 1e9)
#line2, = ax[1].semilogy(sj_x, GetDistribution(0), 'o-')
#ax[1].set_ylim(10^9, 0)

ax.set_xlabel("j")
#ax[1].set_xlabel("s_j")

ax.set_ylabel("N_j")
#ax[1].set_ylabel("N_sj")

#x = np.arange(0, 2*np.pi, 0.01)        # x-array
#line, = ax.plot(x, np.sin(x))

vel_hist = []
avg_hist = []


plt.legend(frameon=False, loc = 2)

tex1 = ax.text(-10,1, '<j>', style='italic', weight=100.0,size=20.0)

for t in range(0, myWF.curTime+1):
    vel_hist.append(GetVelocityJ(t) * 1e4)
    avg_hist.append(GetAvgJ(t))

def animate(i):
    yDat = GetDistribution(i)
    line.set_ydata(yDat)  # update the data
    line2.set_xdata([GetAvgJ(i)]*2)
    line3.set_xdata(avg_hist[0:i])
    line3.set_ydata(vel_hist[0:i])
    tex1.set_x(GetAvgJ(i))
    tex1.set_y(max(vel_hist[i], 3))
    return line,line2,line3,tex1,

#Init only required for blitting to give a clean slate.
def init():
    line.set_ydata(np.ma.array(x, mask=True))
    line2.set_xdata(np.ma.array([0.0], mask=True))
    line3.set_xdata(np.ma.array([GetVelocityJ(i)]*2, mask=True))
    return line,line2,line3,tex1,

ani = animation.FuncAnimation(fig, animate, np.append([0]*100,np.append(np.arange(0, myWF.curTime,3),[myWF.curTime]*100)), init_func=init,
    interval=20, blit=True)
plt.show()