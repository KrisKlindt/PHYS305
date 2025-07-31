# -*- coding: utf-8 -*-
"""
Created on Sun Nov 17 14:02:33 2024

@author: krist
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter

def AdvectionCD(L,N,TotalTime,dt,Skip):

    Nsteps = math.ceil(TotalTime/dt/Skip) # Number of time steps
    
    # define grid
    x,dx = np.linspace(0,L,N,retstep='True')
    Time = [ Skip*dt*i for i in range(Nsteps) ]
    # initialize concentration matrix
    C = np.zeros((N,Nsteps))
    
    # define velocity v at each spatial position
    v = np.sin(x)
    
    # define cos(x) array
    cos = np.cos(x)
    
    # define initial condition
    C[:,0] = 2*np.exp(-(x-(2*np.pi))**2)
        
    # initialize first derivative array
    dC = np.zeros(N)
    
    # initialize second derivative array
    dCd2 = np.zeros(N)


     #setup movie stuff
    metadata = dict(title = 'Q1', artist = 'Matplotlib')
    writer = FFMpegWriter(fps=15, metadata=metadata)
    fig1 = plt.figure()
    l, = plt.plot([], [],'b-')
    plt.xlim(0,L)
    plt.ylim(0,4)

    # begin time stepping
    with writer.saving(fig1,"Q1.mp4",100):
        for i in range(1,Nsteps):
            C[:,i] = C[:,i-1]
            for j in range(Skip):
                # define second derivative
                # impose periodic boundary condition (C[0] = C[N-1])
                dCd2[0] = (C[1,i] - 2*C[0,i] + C[N-2,i])/(dx**2)
                dCd2[1:N-1] = (C[2:N,i] - 2*C[1:N-1,i] + C[:N-2,i])/(dx**2)
                dCd2[N-1] = (C[1,i] - 2*C[N-1,i] + C[N-2,i])/(dx**2)
            
                # compute first derivative
                for a in range(N):
                    # impose periodic boundary condition (C[0] = C[N-1])
                    # v = 0 at boundaries (sin(0) = sin(4pi) = 0), so use forward derivative
                    if a == 0:
                        dC[0] = (C[0,i] - C[N-2,i])/dx
                    elif a == N-1:
                        dC[a] = (C[N-1,i] - C[N-2,i])/dx
                    else:    
                        if v[a] >= 0:
                            dC[a] = (C[a,i] - C[a-1,i])/dx
                        else:
                            dC[a] = (C[a+1,i] - C[a,i])/dx
            
                # time step concentration
                C[:,i] = C[:,i] + dt*1*dCd2 - dt*v*dC - cos*dt*C[:,i]
            
            l.set_data(x,C[:,i])
            writer.grab_frame()
            plt.pause(0.01)
            
    cmap = plt.get_cmap('inferno')
    fig2,ax2 = plt.subplots()
    im = ax2.pcolormesh(Time,x,C,shading='flat',cmap=cmap)
    plt.show()
    return Time,C

T,C = AdvectionCD(4*np.pi,50,10,0.01,5)