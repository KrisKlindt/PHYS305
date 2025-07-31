# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 13:56:34 2024

@author: krist
"""

import math
import numpy as np
from matplotlib.pyplot import plot

def EulerMethod(No):
    # set step sizes, change from 0.01-3
    dt = 1.5
    
    # set Num of steps so that it ends at 30
    Num = math.floor(30/dt)
    
    # intialize array
    N = np.array([0 for i in range(Num)], 'float')
    time = np.array([0 + dt*i for i in range(Num)])
    
    # set first value of array to given initial value
    N[0] = No
    
    for i in range(1,Num):
        N[i] = N[i-1] + dt*(N[i-1] - (N[i-1])**2)
        
    plot(time,N)
    
    return N,time
        
A, B = EulerMethod(0.1)