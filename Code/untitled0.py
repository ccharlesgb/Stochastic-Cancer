
import matplotlib.pyplot as plt
import random
import math

for hid in range(1,4):
    Zn = []
    
    Z_Count = 2000
    Batch_Count = hid * 1000
    
    mean = 4.0
    
    
    for i in range(0,Z_Count):
        
        Batch_Sum = 0.0
        for j in range(0,Batch_Count):
            sample = mean * math.log(1.0/random.random())
            Batch_Sum += sample
        
        Zn.append(Batch_Sum / Batch_Count)
        
    plt.hist(Zn, 30)
    plt.show()