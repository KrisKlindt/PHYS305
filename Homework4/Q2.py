# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 17:23:34 2024

@author: krist
"""
import numpy as np
import matplotlib.pyplot as plt

def RandomWalk():
    Steps = 1000
    Stepsize = 0.1
    
    position = np.zeros((Steps,2), "float")
    Disp = np.zeros((Steps,1), "float")
    SquareDisp = np.zeros((Steps,1), "float")
    AvgDisp = 0
    AvgSquareDisp = 0
    
    StepArray = np.array([0 + i for i in range(Steps)])
    
    for i in range(1,Steps):
        # Random angle for each step
        theta = np.random.uniform(0, 2 * np.pi)
        # Update position based on random direction
        position[i,0] = position[i-1,0] + Stepsize * np.cos(theta)
        position[i,1] = position[i-1,1] + Stepsize * np.sin(theta)
        
        Disp[i,0] = np.sqrt(position[i,0]**2 + position[i,1]**2)
        AvgDisp += Disp[i,0]
        AvgSquareDisp += (position[i,0]**2 + position[i,1]**2)
        SquareDisp[i,0] = AvgSquareDisp/i
        
    AvgDisp = AvgDisp/Steps
    AvgSquareDisp = AvgSquareDisp/Steps

    print(AvgDisp)
    print(AvgSquareDisp)
    
    plt.plot(StepArray, SquareDisp[:,0], label = "Square Displacement vs Steps")
    plt.legend()
    
RandomWalk()