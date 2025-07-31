# -*- coding: utf-8 -*-
"""
Created on Sun Sep 15 21:10:07 2024

@author: krist
"""

import math
import numpy as np
from matplotlib.pyplot import plot
    
def CalculateM(m,c,w):
    # Is used to calculate dm/dt at given time using 4th order Runge-Kutta
    dt = 0.002
    
    k1 = dt *(c*(1-m) - 10 * w * m)
    
    k2 = dt * (c*(1-(m+(k1/2)) - 10 * w * (m+(k1/2))))
    
    k3 = dt * (c*(1-(m+(k2/2)) - 10 * w * (m+(k2/2))))
    
    k4 = dt * (c*(1-(m + k3) - 10 * w * (m+k3)))
    
    return k1,k2,k3,k4
    
def CalculateC(m,c,t):
    # Is used to calculate dc/dt at given time using 4th order Runge-Kutta
    dt = 0.002
    
    S = 0
    if t >= 3 and t <= 3.2:
        S = -2
    else:
        S = 0
        
    k1 = dt * (5*m*(1-c) - 1.25*c + S)
    
    k2 = dt * (5*m*(1-(c+(k1/2))) - 1.25*(c+(k1/2)) + S)
    
    k3 = dt * (5*m*(1-(c+(k2/2))) - 1.25*(c+(k2/2)) + S)
    
    k4 = dt * (5*m*(1-(c+(k3))) - 1.25*(c+(k3)) + S)

    return k1,k2,k3,k4
    
def CalculateW(w,m):
    # Is used to calculate dW/dt at given time using 4th order Runge-Kutta
    dt = 0.02
    
    k1 = dt * (0.1 * (1-w) - 4*m*w)
    
    k2 = dt * (0.1 * (1-(w+(k1/2))) - 4*m*(w+(k1/2)))
    
    k3 = dt * (0.1 * (1-(w+(k2/2))) - 4*m*(w+(k2/2)))
    
    k4 = dt * (0.1 * (1-(w+(k3))) - 4*m*(w+(k3)))
    
    return k1,k2,k3,k4

def RungeKutta(Mo, Co, Wo):
    # set step size
    dt = 0.002
    
    # set Num of steps so that it ends at 30
    Num = math.floor(30/dt)
    
    # intialize arrays for values of x, y, and time
    M = np.array([0 for i in range(Num)], 'float')
    C = np.array([0 for i in range(Num)], 'float')
    W = np.array([0 for i in range(Num)], 'float')
    
    time = np.array([0 + dt*i for i in range(Num)])
    
    # set first value of array to given initial value
    M[0] = Mo
    C[0] = Co
    W[0] = Wo
    
    for i in range(1,Num):
        Mk1, Mk2, Mk3, Mk4 = CalculateM(M[i-1], C[i-1], W[i-1])
        M[i] = M[i-1] + (1/6) * (Mk1 + 2*Mk2 + 2*Mk3 + Mk4)
        
        Ck1, Ck2, Ck3, Ck4 = CalculateC(M[i-1], C[i-1], time[i-1])
        C[i] = C[i-1] + (1/6) * (Ck1 + 2*Ck2 + 2*Ck3 + Ck4)
   
        Wk1, Wk2, Wk3, Wk4 = CalculateW(W[i-1], M[i-1])
        W[i] = W[i-1] + (1/6) * (Wk1 + 2*Wk2 + 2*Wk3 + Wk4)
   
    plot(time,M,time,C,time,W)
    
    return M,C,W
        
Marray, Carray, Warray = RungeKutta(0.0114, 0.0090, 0.9374)