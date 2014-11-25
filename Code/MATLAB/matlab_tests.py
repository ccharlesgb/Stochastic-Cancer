# -*- coding: utf-8 -*-

import scipy.io as io

def SaveRunHistory(file_name, xData, yData, xLabel = "Time", yLabel = "Population"):
    

def Save(file_name, xData, yData, xLabel = "x", yLabel = "y", sim = 0, otherDict = 0):
    data = dict()
    data[xLabel] = xData
    data[yLabel] = yData
    
    if sim != 0:
        data["sim_r0"] = sim.r0   
        data["sim_r1"] = sim.r1
        data["sim_r2"] = sim.r2
        data["sim_u1"] = sim.u1
        data["sim_u2"] = sim.u2
    
    io.savemat(file_name, data)