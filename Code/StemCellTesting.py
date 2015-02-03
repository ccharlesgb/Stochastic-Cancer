# -*- coding: utf-8 -*-
"""
Created on Tue Feb 03 14:15:29 2015

@author: Connor
"""

import Gillespie
import StemCellParam

myGill = Gillespie.Gillespie()

myParam = StemCellParam.StemCellParam()

myGill.Hook(myParam) #Hook into sim parameters

myGill.Simulate()

print ("DONE")
