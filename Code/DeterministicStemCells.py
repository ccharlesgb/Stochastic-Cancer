# -*- coding: utf-8 -*-
"""
Created on Tue Feb 03 15:58:20 2015

@author: Connor
"""

import StemCellParam
import matplotlib.pyplot as plt
import MatTools

myParam = StemCellParam.StemCellParam()

myParam.in0 = 2e2
myParam.in1 = 1e10
myParam.im0 = 0
myParam.im1 = 0

myParam.rn = 0.005
myParam.rm = 0.0115

myParam.cn = 0.75e-2
myParam.cm = 0.38e-2

myParam.dn0 = 0.002
myParam.dn1 = 0.213
myParam.dm0 = 0.002
myParam.dm1 = 0.213

myParam.an = 1.065e3
myParam.am = 1.065e3

def Getn0Dot(param):
    return param.n0 * (param.rn * param.GetPhiNormal() - param.dn0)
    
def Getn1Dot(param):
    return param.n0 * param.an - param.n1 * param.dn1
    
def Getm0Dot(param):
    return param.m0 * (param.rm * param.GetPhiCancer() - param.dm0)
    
def Getm1Dot(param):
    return param.m0 * param.am - param.m1 * param.dm1

deltaT = 1.0
curT = 0.0
maxT = 6000.0

myHist = StemCellParam.StemCellHist()

myParam.Reset()
while curT < maxT:
    curT += deltaT
    n0Dot = Getn0Dot(myParam)
    n1Dot = Getn1Dot(myParam)
    m0Dot = Getm0Dot(myParam)
    m1Dot = Getm1Dot(myParam)
    
    myParam.n0 = myParam.n0 + n0Dot * deltaT
    myParam.n1 = myParam.n1 + n1Dot * deltaT
    myParam.m0 = myParam.m0 + m0Dot * deltaT
    myParam.m1 = myParam.m1 + m1Dot * deltaT
    
    myHist.RecordFrame(curT, myParam)
    
plt.subplot(211)
plt.title("Stem Cells")
plt.plot(myHist.tHist, myHist.n0Hist)
plt.plot(myHist.tHist, myHist.m0Hist, '--')
plt.yscale('log')
plt.ylim(10e0, 10e5)

plt.subplot(212)
plt.title("Differentiated Cells")
plt.plot(myHist.tHist, myHist.n1Hist)
plt.plot(myHist.tHist, myHist.m1Hist, '--')
plt.yscale('log')
plt.ylim(10e0, 10e13)

MatTools.SaveRunHistory()
