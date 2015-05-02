# -*- coding: utf-8 -*-

import scipy.io as io

import os
import time
import matplotlib.pyplot as plt
import inspect

def SaveDict2(data, **kwargs):
    dir_name = inspect.stack()[1][1]; #Which file called us?
    k = dir_name.rfind("/") #Find last appearance of /
    dir_name = dir_name[k+1:-3] #Remove .py and the directory
        
    direct = "MATLAB/Data/S2/" + dir_name
    if not os.path.exists(direct):
        os.makedirs(direct)
    
    file_name = dir_name
    prefix = ""
    for key in kwargs:    
        if key == "prefix":
            prefix = kwargs[key] + "_"
            continue
        file_name = file_name + "_{0}_{1}".format(key, kwargs[key])
    file_name = prefix + file_name
    file_name = "{0}/{1}".format(direct,file_name)
    print("Saving File: {0}".format(file_name))
    
    isFile = os.path.isfile("{0}.mat".format(file_name))
    if isFile == 1:
        print("File already exists overwriting.")
    
    io.savemat(file_name, data)

def SaveDict(file_name, data):
    print("SaveDict is old, use SaveDict2")
    date = time.strftime("%d-%m")
    direct = "MATLAB/Data/S2/" + date
    if not os.path.exists(direct):
        os.makedirs(direct)    
    
    file_name = direct + "/" + file_name   
    print("Saving File: {0}".format(file_name))
    io.savemat(file_name, data)

def ExportFigure(file_name, fig_handle):
    plt.figure(fig_handle[0])
    fig = plt.gcf()
    
    data = dict()
    data["FIG_TITLE"] = fig.get_axes()[0].title
    data["FIG_XLABEL"] = fig.xlabel
    data["FIG_YLABEL"] = fig.ylabel
    
    data["FIG_XMIN"] = fig.xlims[0]
    data["FIX_XMAX"] = fig.xlims[1]
    
    data["FIG_YMIN"] = fig.ylims[0]
    data["FIX_YMAX"] = fig.ylims[1]
    
    lines = fig.get_axes().get_lines()
    
    for i in range(0, len(lines)):
        data["FIG_XDATA_{0}".format(i)] = lines[i].xdata
        data["FIG_YDATA_{0}".format(i)] = lines[i].ydata
        
    SaveDict(file_name, data)

def ColourMap(file_name, xCoords, yCoords, zData, xLabel = "x", yLabel = "y", zLabel = 'z', sim = 0):
    data = dict()
    data[xLabel] = xCoords
    data[yLabel] = yCoords

    data[zLabel] = zData
        
    SaveDict(file_name, data)
    

#Save XY Data from many sims
def SaveXYData(file_name, xData, yData, xLabel = "x", yLabel = "y", yError = [], sim = 0, otherDict = 0):
    data = dict()
    data[xLabel] = xData
    data[yLabel] = yData
    
    if len(yError) > 0:
        data[yLabel + "_ERR"] = yError
    
    if otherDict != 0:
        print("Adding other dict")
        for key, value in otherDict.iteritems():
            print(key)
            data[key] = value    
        
    SaveDict(file_name, data)
    
#For histogram data
def SaveHistogramData(file_name, values, xLabel = "x", yLabel = "Counts", sim = 0, otherDict = 0):
    data = dict()
    data[xLabel] = values
    
    for key, value in otherDict.iteritems():
        data["_" + key] = value      
        
    SaveDict(file_name, data)