# -*- coding: utf-8 -*-
"""
Created on Mon Nov 03 20:14:13 2014

@author: Connor
"""
u2 = 0.5

r1Origin = 1.25
r2Origin = (1-u2)*r1Origin

normGrad = -1.0/(1.0 - u2)

theoryX = []
theoryDiag = []
theoryReac = []
theorySaddle = []

maxr1 = 3.0
minr1 = 0.5

mapSize = 10
for ir1 in range(0,mapSize):
    r1 = (float(ir1)/mapSize)*(maxr1-minr1)+minr1
    r2 = (1-u2)*r1
    theoryX.append(r1)
    theoryDiag.append(r2)
    norm = normGrad * (r1 - r1Origin) + r2Origin
    theoryReac.append(norm)
    
plt.plot(theoryX, theoryDiag, 'k')
plt.plot(theoryX, theoryReac, 'k')
plt.show()