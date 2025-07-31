# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 15:16:56 2024

@author: krist
"""

import numpy as np
import matplotlib.pyplot as plt

def BrownianParticles(N,R,L,e,sigma,TotalTime,dt,Skip):
    # simulates the Brownain dynamics of N neutrally-buoyant
    # particles of radius R moving through water at room
    # temperature and interacting via Lennard-Jones potentials.
    # Periodic boundaries are used with a box size of L.
    # parameters
    
    kBT = 0.004 # thermal energy in pN*um
    eta = 0.001 # viscosity of water
    zeta = 6*np.pi*eta*R # drag coefficient
    c = np.sqrt(2*kBT*zeta/dt) # magnitude of stochastic force
    
    Steps = round(TotalTime/dt/Skip)

    # initialize lists to store positions and time
    x = np.zeros((2,N,Steps))
    Time = [ Skip*dt*i for i in range(Steps) ]
    
    # randomly place particles in box of size L
    x[:,:,0] = L*np.random.rand(2,N)
    
    # initialize matrix to store the particle-particle interactions
    Fi = np.zeros((2,N))

    for i in range(1,Steps):
        x[:,:,i] = x[:,:,i-1]
        
        for j in range(Skip):
            # compute the random force
            xi = c*np.random.randn(2,N)
            
            Dx = np.meshgrid(x[0,:,i])
            Dx = Dx - np.transpose(Dx)
            
            Dy = np.meshgrid(x[1,:,i])
            Dy = Dy - np.transpose(Dy)
            
            Dx[Dx<-L/2] = Dx[Dx<-L/2] + L
            Dx[Dx>L/2] = Dx[Dx>L/2] - L
            
            Dy[Dy<-L/2] = Dy[Dy<-L/2] + L
            Dy[Dy>L/2] = Dy[Dy>L/2] - L
            
            r = np.sqrt(Dx**2 + Dy**2)
            
            for a in range(N):
                r[a,a] = 10
                
            FiMag = (12*e/sigma**2)*( (sigma/r)**(14) - (sigma/r)**8 )
            
            FiMag[FiMag <-0.5] = -0.5
            FiMag[FiMag > 0.5] = 0.5
            
            Fi[0,:] = -np.sum(FiMag*Dx,axis=1)
            Fi[1,:] = -np.sum(FiMag*Dy,axis=1)
            
            x[:,:,i] = x[:,:,i] + (dt/zeta)*(xi + Fi)

            # impose periodic boundary condition
            Mask = x[0,:,i]>L
            x[0,Mask,i] = x[0,Mask,i] - L

            Mask = x[0,:,i]<0
            x[0,Mask,i] = x[0,Mask,i] + L
            
            Mask = x[1,:,i]>L
            x[1,Mask,i] = x[1,Mask,i] - L
            
            Mask = x[1,:,i]<0
            x[1,Mask,i] = x[1,Mask,i] + L
            
    # plots the last time step to see if particles aggregated
    plt.plot(x[0,:,Steps-1], x[1,:,Steps-1], 'ob', markersize = 2)
        
    plt.show()
    
e = np.sqrt(6.84*10**(-11)*0.004)
BrownianParticles(500, 0.01, 1, e, 0.01, 10, 0.001, 10)
#BrownianParticles(N, R, L, e, sigma, TotalTime, dt, Skip)
# e = 5 kBT = 5 * 0.004 = 0.02