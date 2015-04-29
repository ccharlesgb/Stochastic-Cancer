# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 12:14:12 2015

@author: Jonny
"""

import TauLeapParam
import TwoSpeciesAnalytical
import TauSolver
import anfixtime
import TwoSpecies
import wright_fisher

populationSize = 100

#set up parameters
myParam = TauLeapParam.Params(2)
myParam.n0[0] = populationSize
myParam.SetUAll(0.1)
myParam.USE_D = False
myParam.r[0] = 1.0
myParam.r[1] = 1.0

#initialise tau leaping
mySolver = TauSolver.Solver(myParam)
mySolver.CacheX0()


#initialise wright fisher solver
myWF = wright_fisher.wright_fisher(cellTypes)
myHist = TauLeapParam.Hist(cellTypes)

myWF.stopAtAppear = 1

myWF.history = myHist
myWF.timeLimit = 1000000
myWF.useApproxTheta = 0
myWF.params = myParam








myAnalFixtime = TwoSpeciesAnalytical.TwoSpeciesAnalytical(myParam)

fixTimeTest = myAnalFixtime.getFixTime(1)

print("Fix Time New Method is: {0}".format(fixTimeTest))



