# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 12:56:31 2015

@author: Jonny
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 15:59:40 2015

@author: Jonny
"""
import numpy as np
from scipy.optimize import curve_fit

def gauss(x, amp, average, sdev):
    return amp*np.exp(-(x-average)**2/(2.*sdev**2))

class gaussian_fitter:
    def __init__(self):
        #input data        
        
        self.in_x_data = []
        self.in_y_data = []        
        #determines whether all data should be returned, or just the parameters
        self.return_data = 0
        
        #output data
        self.fitted_y_data = []                
        self.var_matrix = []        
        #initial guesses with one array to pass to function
        self.g_mean = 0.0
        self.g_sdev = 1.0
        self.g_amp = 0.0        
        self.p0 = []
        #final values
        self.mean = 0.0
        self.sdev = 0.0
        self.amp = 0.0
        self.p = []
        #hook for the fitter, so that the self argument doesnt confuse it (bless it)        
        self.func = gauss        
        
    def reset(self):
        self.fitted_y_data = []
        self.GetGuess()        
        self.p0 = [self.g_amp, self.g_mean, self.g_sdev]        
        print("p0 is {0}".format(self.p0))
    
    def GetGuess(self):    
        self.g_amp = max(self.in_y_data)
        index = self.in_y_data.index(self.g_amp)        
        self.g_mean=self.in_x_data[index]        
        print("self.amp: {0}, max(data): {1}".format(self.g_amp, max(self.in_y_data)) )
    
    def SetData(self, input_x, input_y):
        self.in_x_data = [0]*len(input_x)
        self.in_y_data = [0]*len(input_y)        
        for i in range(0, len(input_x)):
            self.in_x_data[i] = input_x[i]
            self.in_y_data[i] = input_y[i]

    def FitData(self):         
        self.reset()        
        self.p, self.var_matrix = curve_fit(self.func, self.in_x_data, self.in_y_data, p0=self.p0)
        self.amp,self.mean, self.amp = self.p        
        print(self.p)        
        #if set to return the fitted data, output to array        
        
        if(self.return_data != 0):                
            for i in range(0, len(self.in_x_data) ):
                val = self.func(self.in_x_data[i], *self.p)
                if val < 1:
                    val = 0
                self.fitted_y_data.append(val)
        
        
        
        
        
        
        