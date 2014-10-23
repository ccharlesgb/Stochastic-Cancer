# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import SimTools

pulseOff = 1.0
pulseOn = 2.0
pulseWidth = 0.0
maxPulseWidth = 10.0

pulseWavelength = 10.0

maxTime = 100.0
t = 0.0
step = 0.1

dataT = []
dataP = []
dataA = []

while t < maxTime:
    dataT.append(t)
    #pulseWidth = (float(curPoint) / (dataPointCount - 1)) * maxPulseWidth
    pulseVal = SimTools.PulseWave(t, pulseOn - pulseOff, pulseWidth, pulseWavelength) + pulseOff
    dataP.append(pulseVal)
    avgVal = (pulseOn - pulseOff) * (pulseWidth / pulseWavelength) + pulseOff
    dataA.append(avgVal)
    t += step

plt.plot(dataT, dataP)
plt.plot(dataT, dataA)
plt.show()