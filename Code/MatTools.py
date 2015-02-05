# -*- coding: utf-8 -*-

import scipy.io as io

import os
import time
import matplotlib.pyplot as plt

def SaveDict(file_name, data):
    date = time.strftime("%d-%m")
    direct = "MATLAB/Data/S2/" + date
    if not os.path.exists(direct):
        os.makedirs(direct)    
    
    file_name = direct + "/" + file_name    
    '''
    root = Tkinter.Tk()
    root.withdraw()
    root.focus_force()
    
    #tkFileDialog.Directory(master=direct)
    file_path = tkFileDialog.asksaveasfilename()  
    print("FILE PATH IS {0}".format(file_path))
    '''
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