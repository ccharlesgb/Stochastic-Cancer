# -*- coding: utf-8 -*-
"""
Created on Thu Feb 26 13:35:50 2015

@author: Connor
"""

import math

def SweepParameter(curPoint, pointCount, minParam, maxParam):
    return float(maxParam - minParam) * float(curPoint) / (pointCount - 1.0) + minParam
    
def SweepParamaterLog(curPoint,pointCount,minParam,maxParam):
    progress = float(curPoint) / (pointCount - 1.0)
    a = (float(maxParam)/minParam)
    val = minParam * math.pow(a, progress)
    return val
    
    