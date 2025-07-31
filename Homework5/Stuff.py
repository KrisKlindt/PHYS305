# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 17:13:14 2024

@author: krist
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter

def Diffusion(TotalTime,dt,Skip):
    # simulates the diffusion of particles that are initially
    # localized near the center of a container of size 2
    
    # define parameters
    D = 0.1 # diffusion constant
    Steps = round(TotalTime/dt/Skip)# Number of time steps
    
    # define grid
    N = 100
    # number of grid points
    x,dx = np.linspace(-1,1,N,retstep='True' ,dtype='float')
    Time = [ Skip*dt*i for i in range(Steps) ]

    # initialize array for concentration and set initial condition
    C = np.zeros((N,Steps))
    C[:,0] = np.exp(-x**2/0.1**2)
    
    # initialize second derivative array
    dCd2 = np.zeros(N)
    
    # initialize array for concentration and set initial condition
    metadata = dict(title='Diffusion',artist='Matplotlib')
    writer = FFMpegWriter(fps=15, codec='h264', bitrate=3000, metadata=metadata)
    fig1 = plt.figure()
    data, = plt.plot([],[],'b')
    plt.xlim(x[0],x[-1])
    plt.ylim(0,1)

    with writer.saving(fig1,'Diffusion.mp4',100):
        # begin time stepping
        for i in range(1,Steps):
            C[:,i] = C[:,i-1]
            for j in range(Skip):
                # define second derivative
                dCd2[0] = 2*(C[1,i] - C[0,i])/(dx**2)
                dCd2[1:N-1] = (C[2:N,i] - 2*C[1:N-1,i] + C[:N-2,i])/(dx**2)
                dCd2[N-1] = 2*(C[N-2,i] - C[N-1,i])/(dx**2)
                # time step concentration
                C[:,i] = C[:,i] + dt*D*dCd2
            
            # plot movie frame
            data.set_data(x,C[:,i])
            writer.grab_frame()
            plt.pause(0.02)
        # plot kymograph
        cmap = plt.get_cmap('PiYG')
        fig2,ax2 = plt.subplots()
        im = ax2.pcolormesh(Time,x,C[:N-1,:N-1],shading='flat',cmap=cmap)
        plt.show()
        return C,Time
    
Diffusion(10, 0.001, 10)