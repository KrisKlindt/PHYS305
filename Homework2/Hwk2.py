# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 09:03:02 2024

@author: wolg
"""

import matplotlib.pyplot as plt
import numpy as np

def Problem1(dt):
    
    # define the number of steps to take for a total time of 30
    Steps = round(30/dt) + 1
    
    # initialize vectors
    
    x = np.zeros((Steps,1),'float')
    time = np.array([i*dt for i in range(Steps)],'float')
    
    x[0] = 0.1
    
    for i in range(1,Steps):
        
        x[i] = x[i-1] + dt*(x[i-1]-x[i-1]**2)
        
    plt.plot(time,x)
    
    return time,x

def Problem2(a,b,dt):
    
    # define the number of steps to take for a totaltime of 100
    Steps = round(100/dt) + 1
    
    # initialize vectors
    
    x = np.zeros((Steps,1),'float')
    y = np.zeros((Steps,1),'float')
    time = np.array([i*dt for i in range(Steps)],'float')
    
    x[0] = 1.05*(a+b)
    y[0] = 0.95*b/(a+b)**2
    
    for i in range(1,Steps):
        
        x[i] = ( x[i-1] + dt*(a + x[i-1]**2*y[i-1]) )/(1+dt)
        y[i] = y[i-1] + dt*(b - x[i-1]**2*y[i-1])
        
    plt.plot(time,x,time,y)
    
    return time,x,y

def Problem3(alpha):
    
    # define the number of steps to take for a totaltime of 30
    dt = 0.002
    Steps = round(30/dt) + 1
    
    # initialize vectors
    
    m = np.zeros((Steps,1),'float')
    c = np.zeros((Steps,1),'float')
    w = np.zeros((Steps,1),'float')
    time = np.array([i*dt for i in range(Steps)],'float')
    
    m[0] = 0.0114
    c[0] = 0.0090
    w[0] = 0.9374
    
    S = 0
    
    for i in range(1,Steps):
        
        # take first Runge-Kutta step
        
        # define current values for m, c, and w
        m1 = m[i-1]
        c1 = c[i-1]
        w1 = w[i-1]
        
        k1m = dt*( c1*(1-m1) - 10*w1*m1)
        k1c = dt*( 5*m1*(1-c1) - 1.25*c1 + S )
        k1w = dt*( 0.1*(1-w1) - 4*m1*w1 )
        
        # reevaluate S at half time point
        
        if ( (time[i-1]+dt/2>=3) & (time[i-1]+dt/2 <= 3.2)):
            S = alpha
        else:
            S = 0
        
        # take second Runge-Kutta step
        
        m1 = m[i-1] + k1m/2
        c1 = c[i-1] + k1c/2
        w1 = w[i-1] + k1w/2
        
        k2m = dt*( c1*(1-m1) - 10*w1*m1)
        k2c = dt*( 5*m1*(1-c1) - 1.25*c1 + S )
        k2w = dt*( 0.1*(1-w1) - 4*m1*w1 )       
        
        # take third Runge-Kutta step

        m1 = m[i-1] + k2m/2
        c1 = c[i-1] + k2c/2
        w1 = w[i-1] + k2w/2
        
        k3m = dt*( c1*(1-m1) - 10*w1*m1)
        k3c = dt*( 5*m1*(1-c1) - 1.25*c1 + S )
        k3w = dt*( 0.1*(1-w1) - 4*m1*w1 ) 
        
        # reevaluate S at full time point
        
        if ((time[i-1]+dt>= 3) & (time[i-1]+dt<= 3.2)):
            S = alpha
        else:
            S = 0        
        
        # take fourth Runge-Kutta step
        m1 = m[i-1] + k3m
        c1 = c[i-1] + k3c
        w1 = w[i-1] + k3w
        
        k4m = dt*( c1*(1-m1) - 10*w1*m1)
        k4c = dt*( 5*m1*(1-c1) - 1.25*c1 + S )
        k4w = dt*( 0.1*(1-w1) - 4*m1*w1 ) 
 
        # update m, c, and w
        
        m[i] = m[i-1] + ( k1m + 2*k2m + 2*k3m + k4m )/6
        c[i] = c[i-1] + ( k1c + 2*k2c + 2*k3c + k4c )/6
        w[i] = w[i-1] + ( k1w + 2*k2w + 2*k3w + k4w )/6
        
    plt.plot(time,m,time,c,time,w,linewidth=2)
    plt.xlabel('time',fontsize=14)
    plt.ylabel('concentration',fontsize=14)
    
    return time,m,c,w   