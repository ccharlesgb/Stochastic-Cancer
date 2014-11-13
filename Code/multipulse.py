# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import math
import TwoSpecies

#Class for simulating cancer dynamics with the Gillespie algorithm
class Multipulse:
    def __init__(self):
        
        #define general parameters
        self.drug_strength = 0.0
        self.angle = 0.0
        self.init_r1 = 0.0
        self.init_r2 = 0.0        
        
        #define r1 pulse paramters
        self.freq_r1 = 0.0
        self.r1_width = 0.0 
        self.r1_offset = 0.0
        
        #define r2 pulse parameters       
        self.freq_r2 = 0.0
        self.r2_width = 0.0
        self.r2_offset = 0.0
        
    #helper function to return the contribution to total drug from r1
    def Get_r1_amp(self):
        return self.drug_strength*math.cos(self.angle)
    #helper function to return the contribution to total drug from r2
    def Get_r2_amp(self):
        return self.drug_strength*math.sin(self.angle)
   
   
def multiple_pulse(sim):
    global curTime        
    param=sim.pulseParam        
    #set up the r1 pulse function        
    sim.r1=TwoSpecies.PulseWave(sim.curTime, param.Get_r1_amp(), param.freq_r1, offset=param.r1_offset, width_frac = param.r1_width ) + param.init_r1
    #set up the r2 pulse function
    sim.r2=-TwoSpecies.PulseWave(sim.curTime, param.Get_r2_amp(), param.freq_r2, offset=param.r2_offset, width_frac = param.r2_width ) + param.init_r2
    #print("For angle{0} r1 amp {1} r2 amp {2}".format(param.angle,sim.r1, sim.r2))