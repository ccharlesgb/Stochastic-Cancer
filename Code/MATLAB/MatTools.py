# -*- coding: utf-8 -*-

import scipy.io as io

def SaveRunHistory(file_name, sim, xLabel = "Time", yLabel = "Population"):
    data = dict()
    data[xLabel] = sim.tHist
    data[yLabel + "_n0"] = sim.n0Hist
    data[yLabel + "_n1"] = sim.n1Hist
    data[yLabel + "_n2"] = sim.n2Hist
    
    io.savemat(file_name, data)

def SaveXYData(file_name, xData, yData, xLabel = "x", yLabel = "y", sim = 0, otherDict = 0):
    data = dict()
    data[xLabel] = xData
    data[yLabel] = yData
    
    if sim != 0:
        data["sim_r0"] = sim.r0
        data["sim_r2"] = sim.r2   
        data["sim_r1"] = sim.r1
        data["sim_u1"] = sim.u1
        data["sim_u2"] = sim.u2
    
    for key, value in otherDict.iteritems():
        data["_" + key] = value      
        
    io.savemat(file_name, data)