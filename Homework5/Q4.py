# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 17:46:12 2024

@author: krist
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter

def Q4(TotalTime,dt,Skip):
    
    # define parameters
    D = [0.5, 1, 2, 5] # D constant
    Steps = round(TotalTime/dt/Skip)# Number of time steps
    
    # define grid
    N = 100
    # number of grid points
    x,dx = np.linspace(0,10,N,retstep='True' ,dtype='float')
    Time = [ Skip*dt*i for i in range(Steps) ]

    # initialize array for u and set initial condition
    u = np.zeros((N,Steps))
    u[:,0] = 0.5*(1 - np.tanh(x/0.1))
    u[0,0] = 1   
    
    # initialize second derivative array
    dud2 = np.zeros(N)
    
    t1 = []
    t2 = []
    
    for a in range(len(D)):
        flag1 = 0
        flag2 = 0
        for i in range(1,Steps):
            u[:,i] = u[:,i-1]
            for j in range(Skip):
                # define second derivative
                dud2[0] = 2*(u[1,i] - u[0,i])/(dx**2)
                dud2[1:N-1] = (u[2:N,i] - 2*u[1:N-1,i] + u[:N-2,i])/(dx**2)
                dud2[N-1] = 2*(u[N-2,i] - u[N-1,i])/(dx**2)
                # time step
                u[:,i] = u[:,i] + dt*D[a]*dud2 + (u[:,i-1]*(1-u[:,i-1]))*dt
                if u[49,i] > 0.01 and flag1 != 1: #Checks when x=0.5 has a u value > 0.01
                    t1.append(dt * i * Skip)
                    print(t1)
                    flag1 = 1
                if u[99,i] > 0.01 and flag2 != 1: #Checks when x=1 has a u value > 0.01
                    t2.append(dt*i*Skip)
                    print(t2)
                    flag2 = 1
            
    v = []
    for b in range(len(t1)):
        v.append(0.5/(t2[b]-t1[b])) # v = dx/dt
    
    plt.plot(D, v, label = "Velocity as a function of D")
    plt.legend()
    
Q4(20, 0.001, 10)