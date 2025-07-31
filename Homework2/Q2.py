# -*- coding: utf-8 -*-
"""
Created on Sun Sep 15 19:43:51 2024

@author: krist
"""

import math
import numpy as np
from matplotlib.pyplot import plot

def SemiImplicit(Xo, Yo, A, B):
    # set step size
    dt = 0.1
    
    # set Num of steps so that it ends at 100
    Num = math.floor(100/dt)
    
    # intialize arrays for values of x, y, and time
    X = np.array([0 for i in range(Num)], 'float')
    Y = np.array([0 for i in range(Num)], 'float')
    time = np.array([0 + dt*i for i in range(Num)])
    a = A
    b = B
    
    # set first value of array to given initial value
    X[0] = Xo
    Y[0] = Yo
    
    for i in range(1,Num):
        X[i] = (X[i-1] + dt*(a + X[i-1]**2 * Y[i-1])) / (1 + dt)
        Y[i] = Y[i-1] + dt*(b - X[i-1]**2 * Y[i-1])
        
    plot(time,X,time,Y)
    
    return X,Y
        
Xarray, Yarray = SemiImplicit(0, 0, 1, 2)